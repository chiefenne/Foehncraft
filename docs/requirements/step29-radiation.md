# Step 29: Monte Carlo radiation, version 1 (L15), specification RD0

Specification of 2026-11-24 for step 4 of the recommended order of the
[graphics-algorithms assessment](../assessments/graphics-algorithms-for-modules.md#10-recommended-order)
(the maintainer's "go your recommended order", 2026-11-24). It turns L15 (1)
and the design constraints of D107 into stages that can be built one pull
request at a time. The five questions of section 1 were answered by the
maintainer on 2026-11-24 and are decisions D255 to D259. Coding of each
stage starts on the maintainer's explicit go. Clean-room rules (D201) apply:
methods from the open literature, reference data with a D201 entry before it
appears in a table. Reacting flows are not part of this step; hot gas is given
by its temperature and absorption coefficient.

Sources: the radiation assessment of 2026-11-16 (R1, D107), section 4 of the
graphics-algorithms assessment (null-collision on a majorant grid, stratified
emission, ratio tracking, multiple importance sampling), and the spray
specification of step 28 (the particle layer that radiation shares, D107 (1)).

## 1. Questions and decisions

| # | Question | Decision | Why |
|---|---|---|---|
| Q1 | Lead application, which decides the first gallery case and the validation order | **Answered 2026-11-24 as proposed (D255):** **a box furnace with hot grey gas and cooled tubes** (absorbing-emitting gas, optical thickness from thin to thick by one key), then the heat exchanger of application 203 for surface-to-surface exchange through transparent air | The furnace exercises everything R1 does (gas emission, absorption, walls); exact slab solutions and published cube-enclosure results exist for it. The heat exchanger is the thin-gas benchmark case named in the assessment |
| Q2 | Standalone kernel before coupling | **Answered 2026-11-24 as proposed (D256): yes.** R1a solves radiation on prescribed temperature and absorption fields and is verified against exact solutions; coupling to the energy equation is R1b | The Monte Carlo noise, the ledger closure and the estimators are checked without a flow solve in the loop; R1b then only adds the coupling |
| Q3 | Shared particle layer (D107 (1)) and the order against sprays SP1 | **Answered 2026-11-24 as proposed (D257): the cell walk is shared now, the bundle store later.** R1a adds a walk that reports every cell and its exact path length (`Locator::walk`), built on the tracker that SP1 also uses; the structure-of-arrays bundle store and migration come with SP1 and SP3 and radiation moves onto them | Rays need no equation of motion, only the walk; neither module waits for the other, and there is still one tracker |
| Q4 | Serial first, as sprays (D241) | **Answered 2026-11-24 (D258): multi-rank from R1a**, not the proposed serial start. R1a runs on 1 to N ranks with threads on each rank. A bundle that enters a halo cell is handed to the rank that owns the cell (lengths first, then the bundle records, over `send_recv_f64`), in rounds until the global count of live bundles (`all_sum_i64`) is zero. The number of bundles of an emitter follows from its own power and the exact total power (`sum_exact`), so it does not depend on the partition; the fixed-point tallies are summed as integers. The result is therefore bitwise the same on 1 to 4 ranks. This moves the migration between ranks of D257 forward: R1a builds it on the shared walk, and SP3 uses it for parcels | The maintainer's choice. It builds on MPI W1 to W4 (canonical geometry, halo cells, exact reductions); the standalone kernel needs no flow on ranks. Bundles are keyed by identity (D107 (2)) and the geometry is canonical, so a bundle takes the same path in the same bits on every partition |
| Q5 | GPU from the first pull request, as sprays (D245) | **Answered 2026-11-24 (D259): yes, GPU from R1a**, not the proposed start in R1b. The walk and the radiation kernel are single-source from R1a, compile to PTX in the existing CI job and have a CPU = GPU check (to the tolerance of the exact-sum policy of 2026-11-23, as D245; the integer tallies add no difference of their own). One GPU per run until GPU with ranks comes in version 1.1. Whichever of R1a and sprays SP1 comes first writes the walk on the device; the other reuses it (D257) | The maintainer's choice. Bundles are independent, so the kernel maps directly onto device threads; the device walk is written once for both modules |

## 2. Stages

| Stage | Content | Starts when |
|---|---|---|
| R1a | **Kernel and verification.** Grey absorbing-emitting medium (no scattering) on prescribed T and absorption coefficient per cell; black or grey diffuse walls with a wall temperature, symmetry planes as specular mirrors, open boundaries where bundles escape. Emission: stratified inverse-CDF selection of the emitter (cell or wall face) by emitted power, Owen-scrambled Sobol points (Burley 2020) for position and direction, keyed by (emitter, bundle, call, event) (D107 (2)). Transport: `Locator::walk` with exact path length per cell; absorption weighting with the track-length estimator (each cell takes `w (1 − exp(−κ ℓ))`), the remainder followed to a wall or out; a bundle ends below a weight threshold by Russian roulette. Tallies in 64-bit fixed point (D107 (3)); emitted = absorbed + escaped per call to the fixed-point resolution (D107 (4)). On ranks (D258): hand-over of bundles at halo cells, rounds to a global zero count, bitwise the same tallies on 1 to 4 ranks. On the GPU (D259): single-source walk and kernel, CPU = GPU check. Output: source per cell (W/m³), incident and net flux per wall face, statistical error per cell. Benchmarks (threads, ranks and GPU) and the first gallery images | Q1 to Q5 answered, go |
| R1b | **Coupling to the energy equation.** Source with the emission linearised in T (Fleck & Cummings 1971, D107 (5)), call every N outer iterations with a running average over calls, statistical error as a health number (D69), net wall flux into CHT walls (D160 interface), `[radiation]` case keys, ledger rows of the radiation module (needs 7.0a, D119), on ranks with the energy equation and on the GPU. Gallery case, case file, recipe | R1a, 7.0a, MPI W5 (energy equation on ranks) |
| R1c | Scattering (isotropic, Henyey–Greenstein) and diffuse reflection at grey walls; next-event estimation towards walls with ratio tracking (Novák et al. 2014) | R1b |
| R1d | Reciprocal estimator of the net exchange (Tessé et al. 2002); spike: multiple importance sampling between the forward and the reciprocal estimator (Veach & Guibas 1995) | R1c |
| R1e | Null-collision tracking on a majorant grid (Woodcock et al. 1965; Galtier et al. 2013) as an option; **benchmark against the walk** on the furnace and the heat exchanger (same noise, time per radiation call); the faster one becomes the default | R1a (independent of R1b to R1d) |
| R1f | Weighted sum of grey gases for H₂O/CO₂ mixtures, one grey gas sampled per bundle | R1b |
| R1g | Bundles move onto the structure-of-arrays store of SP1 (D257); the migration between ranks is the one of R1a (D258), which SP3 reuses for parcels | SP1 |

Later, outside version 1: FV-P1 and FV-DOM references (L15 (4)),
full-spectrum k-distributions, the sphere steps of L15 (2), soot and particle
radiation with sprays.

## 3. Validation list (exact or open literature)

R1a:

1. **Slab between black plates**, isothermal grey medium, cold walls, side
   walls as symmetry planes (an infinite slab): wall flux
   `q = σ T⁴ (1 − 2 E₃(τ))` with the exponential integral `E₃` and the optical
   thickness `τ` (Modest & Mazumder 2021, chapter on the one-dimensional
   slab), for τ from 0.01 to 10; the error falls as 1/√N and stays inside
   three standard errors (D90 amendment).
2. **Isothermal enclosure**: walls and gas at one temperature: net source per
   cell and net wall flux zero within three standard errors; the reciprocal
   estimator of R1d makes it exactly zero.
3. **Transparent medium, two black parallel plates**: net exchange against
   the exact view factor of parallel rectangles (Howell's catalogue of
   configuration factors).
4. **Ledger**: emitted = absorbed + escaped to the fixed-point resolution on
   every call; bitwise the same tallies on 1, 2 and 4 threads (D93) and on
   1, 2 and 4 ranks (D258, in the MPI CI job); CPU = GPU (D259).
5. **Cube enclosure with an isothermal grey gas** and cold black walls: wall
   flux distribution against published benchmark solutions (D201 entry first).
6. **Sobol against pseudo-random emission**: error against N on cases 1
   and 3.

R1e: walk against null-collision at equal statistical error (time per
call), on the furnace and on the heat exchanger.

## 4. Benchmark rule

By the maintainer's rule for geometry topics (2026-11-21), every stage
names the graphics method it is measured against: the walk against the
null-collision tracker (R1e), the track-length estimator against the
collision estimator, Sobol against pseudo-random points, and the forward
against the reciprocal and the MIS estimator (R1d). Embree stays a timing
reference only (float32; assessment section 4.3).

## 5. Sizing

R1a: four pull requests (walk with the hand-over between ranks; kernel
with cases 1 to 4 on threads and ranks; walk and kernel on the GPU with the
CPU = GPU check; output, cube case and gallery), about 45 to 60 minutes of
code each, verification runs of minutes.

Status of R1a (coding started on the maintainer's go, 2026-11-24):

1. **Walk and hand-over between ranks: #254.** The walk body
   `foehncraft_kernels::walk` (one source for CPU and GPU; the particle
   tracker of D47 now runs on it, with bitwise its former results) reports
   every cell with the interval of the segment inside it;
   `foehncraft_mesh::handover` hands a walk stopped at a processor face to
   the rank behind it. Tests: the intervals tile the path, the exit point
   is the first hit of the ray cast of D108 part 2 to 10⁻¹², the lengths in
   a hex mesh are the chords; on 2 to 4 ranks (RCB and METIS, tetrahedral,
   polyhedral and mixed meshes) every interval and every end are bitwise
   those of one rank, also with MPI in the CI job (`examples/mpi_walk.rs`,
   1 to 4 ranks, up to 2 311 hand-overs per run).
2. **Kernel with the cases 1 to 4 on threads and ranks: #255.** The bundle
   kernel `foehncraft_kernels::radiation` (scrambled Sobol emission,
   track-length absorption, grey walls, mirrors, open boundaries, roulette,
   integer quanta) and the crate `foehncraft-radiation` (emitters, rounds
   of hand-overs, ledger, batch errors, exact solutions). Case 1 within
   0.14 % of the exact wall flux for τ = 0.01 to 10 and every layer source
   within three standard errors; case 2 with a mean net source of about 2 %
   of the emission; case 3 within one standard error of the configuration
   factor; case 4 closes the ledger to the quantum and gives bitwise the
   same tallies on 1, 2 and 4 threads and on 1 to 4 ranks (also with MPI in
   CI, `examples/mpi_radiation.rs`). The batch-mean error is conservative:
   the bundles of an emitter are stratified, the batches are not.
3. **Walk and kernel on the GPU** (in review). The PTX entry
   `radiation_transport` (`gpu/ptx-kernels/src/radiation.rs`) runs the
   kernel body of the CPU for one bundle per GPU thread; the tallies are
   whole numbers of quanta added with 64-bit integer atomics, the same sum in
   any order (D269, proposed: the one exception of R10, checked by
   `tools/check/kernels.py`, which allows integer adds only in this entry).
   The bundles are emitted on the host and sent in launches of 2²⁰. The
   CUDA backend (`crates/foehncraft-cuda/src/radiation.rs`) registers with
   `foehncraft_radiation::register_gpu`; `Domain::on_gpu` uses it. The data
   path of the GPU (flat records, the entry body, tallies per launch) runs
   on the CPU in the tests and gives bitwise the CPU result. The CPU = GPU
   check on a GPU: `python3 -m tools.bench.radiation_cloud` or
   `tools/notebooks/radiation_colab.ipynb`.
4. Output, cube case and gallery. R1b to R1g: not sized yet.

## 6. References

Burley 2020 (Owen scrambling, JCGT 9(4)); Fleck & Cummings 1971 (J. Comput.
Phys. 8); Galtier et al. 2013 (JQSRT 125); Modest & Mazumder 2021 (Radiative
Heat Transfer, 4th ed.); Howell, Mengüç, Daun & Siegel 2021 (Thermal Radiation
Heat Transfer, 7th ed.); Novák, Selle & Jarosz 2014 (residual ratio tracking,
ACM TOG 33(6)); Tessé et al. 2002 (Int. J. Heat Mass Transfer 45); Veach &
Guibas 1995 (SIGGRAPH); Woodcock et al. 1965 (ANL-7050).
