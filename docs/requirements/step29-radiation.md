# Step 29: Monte Carlo radiation, version 1 (L15), specification RD0

Draft of 2026-11-24 for step 4 of the recommended order of the
[graphics-algorithms assessment](../assessments/graphics-algorithms-for-modules.md#10-recommended-order)
(the maintainer's "go your recommended order", 2026-11-24). It turns L15 (1)
and the design constraints of D107 into stages that can be built one pull
request at a time. The decisions in section 1 are **proposals until the
maintainer answers them**; each gets its D number when answered. Coding of each
stage starts on the maintainer's explicit go. Clean-room rules (D201) apply:
methods from the open literature, reference data with a D201 entry before it
appears in a table. Reacting flows are not part of this step; hot gas is given
by its temperature and absorption coefficient.

Sources: the radiation assessment of 2026-11-16 (R1, D107), section 4 of the
graphics-algorithms assessment (null-collision on a majorant grid, stratified
emission, ratio tracking, multiple importance sampling), and the spray
specification of step 28 (the particle layer that radiation shares, D107 (1)).

## 1. Questions and proposed answers

| # | Question | Proposal | Why |
|---|---|---|---|
| Q1 | Lead application, which decides the first gallery case and the validation order | **Answered 2026-11-24 as proposed (D247):** **a box furnace with hot grey gas and cooled tubes** (absorbing-emitting gas, optical thickness from thin to thick by one key), then the heat exchanger of application 203 for surface-to-surface exchange through transparent air | The furnace exercises everything R1 does (gas emission, absorption, walls); exact slab solutions and published cube-enclosure results exist for it. The heat exchanger is the thin-gas benchmark case named in the assessment |
| Q2 | Standalone kernel before coupling | **Yes.** R1a solves radiation on prescribed temperature and absorption fields and is verified against exact solutions; coupling to the energy equation is R1b | The Monte Carlo noise, the ledger closure and the estimators are checked without a flow solve in the loop; R1b then only adds the coupling |
| Q3 | Shared particle layer (D107 (1)) and the order against sprays SP1 | **The cell walk is shared now, the bundle store later.** R1a adds a walk that reports every cell and its exact path length (`Locator::walk`), built on the tracker that SP1 also uses; the structure-of-arrays bundle store and migration come with SP1 and SP3 and radiation moves onto them | Rays need no equation of motion, only the walk; neither module waits for the other, and there is still one tracker |
| Q4 | Serial first, as sprays (D241) | **Yes.** R1a to R1d run on one process (threads allowed); a radiation case on more than one rank is refused with a diagnostic code; migration between ranks after MPI W3 to W5 | The same rule as sprays; bundles are keyed by identity (D107 (2)), so nothing is rewritten for ranks |
| Q5 | GPU from the first pull request, as sprays (D245) | **No, from R1b.** R1a is CPU; the GPU kernels come once SP1 has put the walk on the device, with the CPU = GPU check | The walk on the GPU is SP1's work; doing it twice wastes it. Bundles are independent, so the CPU code is already in the form the wavefront kernels need |

## 2. Stages

| Stage | Content | Starts when |
|---|---|---|
| R1a | **Kernel and verification.** Grey absorbing-emitting medium (no scattering) on prescribed T and absorption coefficient per cell; black or grey diffuse walls with a wall temperature, symmetry planes as specular mirrors, open boundaries where bundles escape. Emission: stratified inverse-CDF selection of the emitter (cell or wall face) by emitted power, Owen-scrambled Sobol points (Burley 2020) for position and direction, keyed by (emitter, bundle, call, event) (D107 (2)). Transport: `Locator::walk` with exact path length per cell; absorption weighting with the track-length estimator (each cell takes `w (1 − exp(−κ ℓ))`), the remainder followed to a wall or out; a bundle ends below a weight threshold by Russian roulette. Tallies in 64-bit fixed point (D107 (3)); emitted = absorbed + escaped per call to the fixed-point resolution (D107 (4)). Output: source per cell (W/m³), incident and net flux per wall face, statistical error per cell. Benchmarks and the first gallery images | Q1 to Q5 answered, go |
| R1b | **Coupling to the energy equation.** Source with the emission linearised in T (Fleck & Cummings 1971, D107 (5)), call every N outer iterations with a running average over calls, statistical error as a health number (D69), net wall flux into CHT walls (D160 interface), `[radiation]` case keys, ledger rows of the radiation module (needs 7.0a, D119). GPU kernels with the CPU = GPU check. Gallery case, case file, recipe | R1a, 7.0a |
| R1c | Scattering (isotropic, Henyey–Greenstein) and diffuse reflection at grey walls; next-event estimation towards walls with ratio tracking (Novák et al. 2014) | R1b |
| R1d | Reciprocal estimator of the net exchange (Tessé et al. 2002); spike: multiple importance sampling between the forward and the reciprocal estimator (Veach & Guibas 1995) | R1c |
| R1e | Null-collision tracking on a majorant grid (Woodcock et al. 1965; Galtier et al. 2013) as an option; **benchmark against the walk** on the furnace and the heat exchanger (same noise, time per radiation call); the faster one becomes the default | R1a (independent of R1b to R1d) |
| R1f | Weighted sum of grey gases for H₂O/CO₂ mixtures, one grey gas sampled per bundle | R1b |
| R1g | Migration between ranks (with SP3) and the multi-rank check (bitwise the same on 1 to 4 ranks) | MPI W3 to W5, SP3 |

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
   every call; bitwise the same tallies on 1, 2 and 4 threads (D93).
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

R1a: two pull requests (walk and kernel with cases 1 to 4; output, cube
case and gallery), about 45 to 60 minutes of code each, verification runs of
minutes. R1b to R1g: not sized yet.

## 6. References

Burley 2020 (Owen scrambling, JCGT 9(4)); Fleck & Cummings 1971 (J. Comput.
Phys. 8); Galtier et al. 2013 (JQSRT 125); Modest & Mazumder 2021 (Radiative
Heat Transfer, 4th ed.); Howell, Mengüç, Daun & Siegel 2021 (Thermal Radiation
Heat Transfer, 7th ed.); Novák, Selle & Jarosz 2014 (residual ratio tracking,
ACM TOG 33(6)); Tessé et al. 2002 (Int. J. Heat Mass Transfer 45); Veach &
Guibas 1995 (SIGGRAPH); Woodcock et al. 1965 (ANL-7050).
