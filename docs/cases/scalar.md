# Scalar cases (`examples/cases/scalar/`)

Part of the [Case catalogue](../cases.md): steady convection–diffusion of a scalar: manufactured solutions, convection schemes; scalar keys.

On this page:

- [Scalar case-file keys](#scalar-case-file-keys)

Stage 2: steady convection–diffusion of a scalar, ∇·(ρuφ) − ∇·(Γ∇φ) = s,
with the operators of stage 2 ([explanation](../verification/stage2.md)). Run
with `foehn scalar FILE`; command-line options override the file
(`foehn scalar --help`). Each run writes `scalar.md` and `scalar.json`
(per solve: iterations, residual, error norms, φ range, time per phase, ns
per cell and iteration, comparison with the first backend), `fields.vtkhdf`
(φ, limiter, error), `convergence.svg` for a sweep and `profile.svg` for a
sampling line. The run fails when a check of the case fails (observed
order, bounds, convergence, backends differing by more than 1e-10). CI runs
every case except `*_gpu.toml` (R5).

| Case file | Mesh | Shows | Where |
|---|---|---|---|
| `scalar/mms_poly_box.toml` | generated poly 3D, 729 to 35 937 cells (3 levels) | manufactured solution, order of accuracy (1.8–2.2 required) | [gallery 8](../gallery/scalar-pressure.md#scalar-mms), [recipe 12](../recipes/scalar-pressure.md#scalar-mms) |
| `scalar/mms_hex_box.toml` | generated hex 3D, 512 to 32 768 cells | the same on perturbed (warped) hexahedra | [recipe 12](../recipes/scalar-pressure.md#scalar-mms) |
| `scalar/mms_tet_box.toml` | generated tet 3D, 3 072 to 196 608 cells | the same on tetrahedra | [recipe 12](../recipes/scalar-pressure.md#scalar-mms) |
| `scalar/mms_mixed_box.toml` | generated mixed 3D, 1 344 to 95 232 cells | the same on hexahedra, prisms, tetrahedra and polyhedra | [recipe 12](../recipes/scalar-pressure.md#scalar-mms) |
| `scalar/mms_poly_2d.toml` | generated poly 2D, 289 to 16 641 cells (4 levels) | planar solution, zero-gradient (symmetry) z boundaries | [recipe 12](../recipes/scalar-pressure.md#scalar-mms) |
| `scalar/oblique_step_2d.toml` | generated poly 2D, 2 401 cells | upwind, linear upwind, both limiters and the bounded central scheme (stage 5, D25) as `[[variants]]`; bounds check | [gallery 9](../gallery/scalar-pressure.md#oblique-step), [recipe 13](../recipes/scalar-pressure.md#scalar-schemes) |
| `scalar/mms_poly_gpu.toml` | generated poly 3D, 912 673 cells | CPU against two GPU runs: GPU equals CPU and GPU run to run | [recipe 14](../recipes/gpu.md#scalar-gpu), `tools/notebooks/stage2_colab.ipynb` |
| `scalar/ranks_oblique_step_2d.toml` | generated poly 2D, 4 225 cells | the oblique step with linear upwind and the Venkatakrishnan limiter on 2, 3 and 5 ranks, bitwise the one-rank solution (MPI step W3) | [gallery 230](../gallery/meshes.md#scalar-ranks), [recipe 230](../recipes/meshes.md#scalar-ranks) |
| `scalar/poisson_amg_poly_box.toml` | generated poly 3D, 729 to 35 937 cells | Poisson manufactured solution with CG + AMG (stage 3) against BiCGStab + Jacobi as `[[variants]]`; order of accuracy | [gallery 11](../gallery/scalar-pressure.md#poisson-amg-mms), [recipe 15](../recipes/scalar-pressure.md#poisson-study) |

## Scalar case-file keys

`[mesh]` and `[mesh.generate]` as in the [mesh cases](meshes.md#mesh-case-file-keys).
A boundary condition is required for every boundary.

| Key | Default | Meaning |
|---|---|---|
| `problem.density` | 1.0 | ρ (kg/m³) |
| `problem.diffusivity` | 0.05 | Γ (kg/(m s)) |
| `problem.velocity` | `[1.0, 0.6, 0.3]` | u (m/s), constant |
| `problem.source` | 0 | constant volumetric source |
| `problem.manufactured` | none | `{}` for the default manufactured solution, or `{ offset, amplitude, wavenumber, phase }`; sets the source, the `manufactured` boundary values and the error norms |
| `problem.boundary.<boundary>.type` | – (required) | `fixed-value`, `fixed-gradient` (outward normal gradient) or `zero-gradient` (also symmetry and outflow) |
| `problem.boundary.<boundary>.profile` | `"constant"` | `constant` (`value`), `step` (`axis`, `at`, `below`, `above`), `linear` (`value`, `gradient`), `manufactured` |
| `solver.scheme` | `"linear-upwind"` | `upwind`, `linear-upwind` or `central` (bounded central, D25) |
| `solver.central_blending` | 1 | with `central`: share c of the bounded central correction; the rest, 1 − c, is linear upwind |
| `solver.central_bounding` | 0.1 | with `central`: β of the normalised-variable bound (central where φ̃_U ≥ β, blended towards upwind below, upwind outside (0, 1)); 0 for plain central |
| `solver.limiter` | `"none"` | `none`, `venkatakrishnan`, `venkatakrishnan-wang`, `barth-jespersen` (ψ = min(1, Δ₊/Δ₋), strictly bounded) |
| `solver.limiter_k` | 5.0 | K: ε² = (K h)³ (Venkatakrishnan) or ε = K (φ_max − φ_min) (Venkatakrishnan–Wang); 0 for strictly bounded |
| `solver.nonorth_correction` | true | non-orthogonal diffusion correction |
| `solver.nonorth_limit` | −1 | limit of that correction as a multiple of the orthogonal part (negative: none) |
| `solver.relaxation` | 1.0 | implicit under-relaxation factor |
| `solver.max_outer` | 200 | outer (deferred-correction) iterations |
| `solver.tolerance` | 1e-10 | outer tolerance on ‖b − Aφ‖/‖b‖ |
| `solver.linear_max_iter` | 200 | linear iterations per outer iteration |
| `solver.linear_tolerance` | 1e-3 | relative tolerance of each linear solve |
| `solver.linear_solver` | `"bicgstab"` | `bicgstab` (any matrix) or `cg` (flexible CG, symmetric positive definite matrices: no convection) |
| `solver.preconditioner` | `"jacobi"` | `jacobi` or `amg` (one V-cycle of the smoothed-aggregation multigrid, stage 3) |
| `solver.amg` | defaults | multigrid settings, the keys of [`linear.amg`](poisson.md#poisson-case-file-keys) |
| `run.backends` | `["cpu"]` | `cpu`, `cpu-serial`, `cuda`; later backends are compared with the first |
| `run.threads` | all cores | CPU threads |
| `run.rust_ptx` | `target/nvptx64-nvidia-cuda/release/foehncraft_ptx_kernels.ptx` | PTX module of the GPU backend |
| `run.ranks` | `[]` | also solve every variant on these numbers of ranks, as threads of this process (RCB partition, MPI step W3); every solve must be bitwise the first backend's, else the run fails; Jacobi preconditioner only until W4 |
| `verification.levels` | none | refinement levels of `[mesh.generate]`; two or more: order-of-accuracy sweep |
| `verification.expected_order` | none | `[lo, hi]`: fail if the order on the finest pair is outside |
| `verification.bounds` | none | `[lo, hi]`: fail if φ leaves this range (by more than 1e-6) |
| `[[variants]]` | none | `name`, `solver = { ... }`, optional `bounds`: the case once per variant |
| `output.directory` | `output/<case name>` | output directory |
| `output.vtkhdf` | true | write `fields.vtkhdf` |
| `output.profile` | none | `{ from, to, points }`: sample φ along a line (`profile.csv`, `profile.svg`) |
