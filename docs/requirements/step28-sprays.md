# Step 28: Lagrangian sprays, version 1 (L12), specification SP0

Decisions of the maintainer, 2026-11-24, from the assessment
[spray-module-prerequisites](../assessments/spray-module-prerequisites.md)
(PR #228). Specification only; coding of each part starts on the maintainer's
explicit go. Clean-room rules (D201) apply: models from the open literature,
constants as case keys with the literature value as default, no numbers copied
from commercial documentation. Reacting flows are not part of this step.
The validation matrix, the statistics and datum reports and the MurmelCraft
viewer (D261 to D268) are in
[step28-spray-validation-statistics](step28-spray-validation-statistics.md).

## 1. Decisions

| ID | Decision |
|---|---|
| D240 | **Lead application: water-spray cooling of a hot gas stream.** It decides the validation list and the first gallery case: injectors spraying water into a hot air duct, gas cooled by sensible heat first (SP2), then by evaporation (SP5). Urea-water for SCR, fuel sprays, painting and scrubbers follow with their own cases. |
| D241 | **Serial release first.** SP0 to SP2, SP4 and SP5 run on one process (threads allowed). A spray case started on more than one rank is refused with a diagnostic code. Parcel migration between ranks (SP3) is planned after MPI W3 to W5 are merged. All parcel kernels are pure functions of (parcel, gas state), so nothing is rewritten for SP3. |
| D242 | **Droplet size input.** Three forms as case keys: constant diameter, Rosin-Rammler (Rosin and Rammler 1933) and log-normal; a table form is not in v1. SP1 ships constant diameter, SP2 adds the two distributions. Sizes are drawn from the counter-based generator (stream `SPRAY_INJECTION`, D121). |
| D243 | **Turbulent dispersion: eddy-lifetime model first** (Gosman and Ioannides 1983), using k and epsilon (or omega) at the parcel, stream `SPRAY_DISPERSION`. **Answered by the maintainer 2026-11-24: both models as values of one case key** (`dispersion = "eddy-lifetime"` or `"random-walk"`, same stream): eddy lifetime ships in SP2, the continuous random walk follows behind the same key. |
| D244 | **The spray ledger pair waits for the ledger as data (D119, stage 7 part 7.0a).** SP2 (two-way coupling) starts only when 7.0a is on `main`; the gas-parcel exchange (momentum, mass, enthalpy) is then booked as ledger rows with opposite signs on both sides. Consequence: 7.0a becomes a prerequisite of SP2, not of SP1. **Done 2026-11-24 (PR 238):** `Rates::book`, storage and closure for every account, thermal account, ledger in `result.json`; a closed-box test with a stand-in module closes mass, momentum and heat to round-off. |
| D245 | **GPU kernels from SP1.** Parcel kernels are single-source from the first pull request, compile to PTX in the existing CI job, and have a CPU = GPU check (to the tolerance of the exact-sum policy of 2026-11-23). SP1 is therefore larger than planned; the separate stage SP9 reduces to multi-GPU and performance work. |
| D246 | **Injector keys in an included file** (D114). The case file names the spray file (`[spray] include = "spray.toml"`); injectors (`[[injector]]`: position, direction or cone, mass flow over time, size form, liquid, temperature) and size or rate tables live there. The resolved case (`case.resolved.toml`) holds the merged keys, so a run stays reproducible and agent-readable. |

## 2. Stages after the decisions

| Stage | Content | Starts when |
|---|---|---|
| SP1 | Particle layer, one rank, one-way: SoA parcels, `[[injector]]`, tracker on `Locator`, exponential sub-stepping, Schiller-Naumann drag, gravity, stick and rebound walls, parcel block in the VTKHDF output and restart, **GPU kernels with CPU = GPU check**. Gallery example, case, recipe. | go |
| 7.0a | Ledger as data (D119), in stage 7 | before SP2 |
| SP1b | Datum object for gas and spray, liquid flow through datums, droplet statistics, structured output (D262 to D267), checks V28.12 to V28.15; details in [step28-spray-validation-statistics](step28-spray-validation-statistics.md) | SP1 merged and go |
| SP2 | Two-way momentum and sensible heat, size distributions, eddy-lifetime dispersion, ledger rows for the gas-parcel pair. The coupling follows section 5 (D247 to D252, D254): frozen gas state, sorted deposit, implicit relaxation form, closing step, passes inside the time step | SP1 and 7.0a |
| SP4 | Minimal vapour species on the low-Mach density, liquid property set for water first | SP2 |
| SP5 | Evaporation (Ranz-Marshall, Abramzon-Sirignano), closed-box mass and energy ledger test, the spray-cooling case. Vapour and enthalpy exchange in the relaxation form of D250, parcel integrators of D251, saturation limiter of D253, checks V28.1 to V28.8 | SP4 |
| SP2s | Steady spray mode (D260): parcel streams through the steady gas field, relaxed sources, checks V28.9 to V28.11. Momentum and sensible heat first; evaporation in the steady mode when SP5 is merged. **Implemented 2026-11-25 (PR #292)** for momentum and sensible heat; see section 6 | SP2 merged and V28.1 to V28.4 passed |
| SP6 | Breakup (KH-ACT and KH-RT with blob injection, TAB distortion with dynamic drag, rate tables, mirrored eddies; D300 to D308, section 7), then collision | owner request 2026-11-25 (ECN Spray G) |
| SP3 | Parcel migration between ranks | MPI W3 to W5 merged |
| SP7 | Wall film | 7t2 surface mesh |

## 3. Validation list (open literature or exact)

Stokes and Schiller-Naumann terminal velocity; exponential velocity
relaxation in a uniform stream; streamline test in a prescribed vortex; bitwise
results on 1 to 4 threads and CPU = GPU; closed-box momentum balance
(gas + parcels + wall to round-off); dispersion in grid turbulence (Snyder and
Lumley 1971, Wells and Stock 1983); d-squared law (Godsave 1953, Spalding
1953); closed-box evaporation with mass and enthalpy ledgers to round-off and
entropy production not negative (D90); the spray-cooling case against the
energy balance of the duct. Every reference data set gets its D201 entry before
it appears in a table. The measured anchors, the own ODE references and the
stage of each check (V28.12 to V28.21) are in
[step28-spray-validation-statistics](step28-spray-validation-statistics.md)
(D261).

## 4. Open points

- D243 (dispersion model) is answered: both models as one key (2026-11-24).
- SP6: aerodynamic breakup of fuel sprays first (owner request 2026-11-25,
  ECN Spray G and Spray A; section 7). Collision, a nozzle-turbulence and
  cavitation term and the fragment-size distributions of slow drops are
  later parts.
- Turbulence modulation by particles is in no decision yet (see the
  validation page, open points).
- Liquids beyond water (ethanol, n-heptane, n-dodecane) need source entries
  in `DATA_SOURCES.md` before use.

## 5. Two-way coupling: independent of parcel order, accurate at practical time steps (D247–D254)

Question of the maintainer, 2026-11-24: is all spray work parallel on CPU and
GPU, and how do we avoid the classic error in which the parcels of one cell
change that cell (vapour, temperature, pressure) while parcels elsewhere,
handled at another moment, see other conditions, so that only small time steps
hide the inaccuracy? This section is the answer and the design. It binds SP2
(momentum and sensible heat) and SP5 (evaporation); SP1 (one-way) keeps the
parts marked "from SP1".

### 5.1 Where the classic error comes from

Three different things are usually mixed up:

1. **Order dependence.** A parcel loop that updates the cell while it runs
   (the next parcel in the cell sees the vapour and the temperature the
   previous one left) gives a result that depends on the order of the parcels,
   and so on the thread count, the partition and the sort. This is an error of
   construction, not of the time step.
2. **Lagged, explicit sources.** The classic particle-source-in-cell method
   (Crowe, Sharma and Stock 1977) freezes the gas, integrates all parcels, and
   hands the summed sources to the next gas solve as fixed numbers. It is
   independent of the order, but each parcel evaporates as if it were alone:
   the sum over a cell can pass saturation, and the gas temperature can pass
   the droplet temperature. As a lagged relaxation it is stable only for
   Δt < 2 τ and free of oscillation only for Δt < τ, where τ is the coupling
   time of the cell (below). This is the reason for the small time steps.
3. **Splitting in time.** Gas and parcels advance one after the other, so the
   exchange belongs to another time level than the transport that carries it
   to the neighbour cells. The error is of first order in Δt and vanishes only
   when the two are solved together.

The coupling times of a cell with droplet number density n, diameter d and
loading Z (liquid mass over gas mass), from the relaxation forms of the
exchange (own estimate from the cited laws, not a reference value):

| Exchange | Coupling time of the gas in the cell | Example: water, d = 20 µm, n = 10¹⁰ m⁻³ (liquid volume fraction 4·10⁻⁵) |
|---|---|---|
| Momentum | τ_p / Z, with τ_p = ρ_l d² / (18 μ) | τ_p about 1.2 ms |
| Vapour (saturation) | τ_sat = (1 − Y_s) / (n π d D Sh) | 10 to 30 ms |
| Heat | τ_T,gas = ρ c_p / (n π d k Nu) | same order as τ_sat |

Near an injector n is 10 to 100 times larger and the times fall below 1 ms,
while a practical flow step of the lead application (D240) is 0.1 to 10 ms.
The design must therefore be right for Δt far above τ, not only below it.

### 5.2 Decisions

| ID | Decision |
|---|---|
| D247 | **Parallelism contract of the spray.** A *parcel pass* is a pure function: (parcel at the old time level, frozen gas state, keyed random numbers) gives (parcel at the new level, its exchange records). No parcel reads what another parcel wrote in the same pass. The pass therefore runs over the parcels in any order: CPU threads and GPU from SP1 (D245), ranks with SP3 (D241). Results (fields, parcels, ledger) are bitwise the same for any parcel storage order and any thread count; for any rank count with SP3; CPU against GPU to the tolerance of D245. Rules that keep it so: (a) sub-steps of a parcel end at cell faces and at the parcel's own accuracy limits (D251), never at a thread chunk, a rank face, an output time or a pass number; a rank face is a cell face like any other; (b) random numbers are keyed by (stream, parcel id, time step, event) and never by the pass, the sub-step counter of another parcel or the thread (D121), so a repeated pass draws the same numbers; (c) models in which parcels see each other (collision, SP6) get a rule of their own before they are built: partners from the sorted list of the cell, both updated from the state before the collision. |
| D248 | **Frozen gas state: all parcels at the same moment.** Every parcel of a pass sees one and the same gas state: the cell values and least-squares gradients of velocity, temperature, vapour fraction, density, the turbulence quantities and the one thermodynamic pressure p₀ of the low-Mach form (D151). It is the state at the new time level in its current iterate (first pass: the state the time step starts from, extrapolated as `solver.extrapolation` does for the explicit terms). The gas is not changed during a pass. Cells and their neighbours are so treated at the same instant; the exchange reaches the neighbours only through the transport equations, which hold all sources at once. p₀ is one number for the whole closed volume, so evaporation in one cell cannot give that cell a thermodynamic pressure of its own; the volume the vapour needs enters the pressure-correction equation as a mass source in the same outer iteration, and p₀(t) follows from the total mass including the vapour. |
| D249 | **Deterministic deposit; one number, two books.** Each sub-step of a parcel writes one record: the cell (global id), the parcel (global id), the sub-step number, the parcel's own change of mass, momentum and enthalpy in that sub-step, and the slopes of D250. The exchange *is* the parcel's change (Crowe, Sharma and Stock 1977), split over the cells by the time spent in each, so nothing is evaluated twice and conservation holds to round-off. Records are sorted by (cell id, parcel id, sub-step) and summed per cell in that order (D117 (2)); no floating-point atomics (D31). The same sorted sums run on the GPU (sort by key, then one sequential sum per cell), so the per-cell sources are bitwise equal on CPU and GPU. The ledger books each pair from the same number with opposite signs (D117, D119, D244). |
| D250 | **Relaxation form with implicit gas sources, and a closing step.** Every exchange is handed to the gas as *value plus slope*: S = E − β (φ − φ*), with φ* the frozen cell value the pass used, E the summed exchange and β ≥ 0 the summed slope. β goes to the diagonal of the gas equation (Patankar 1980), so the gas relaxes towards the state of its parcels and cannot pass it, at any Δt. Momentum: β = N m (1 − e^(−δt/τ_p)) per sub-step, from the exact relaxation of D251 (implicit momentum coupling as in Dukowicz 1980). Heat: the same with the thermal response time, entered on the temperature diagonal of the energy equation (D150). Vapour: the *chord* to the saturation state at the droplet surface, β = E / (Y_s − Y*), not the tangent of ln(1 + B_M): the tangent of that concave law overshoots saturation, the chord ends on it. A parcel that evaporates completely in the step has β = 0 (its mass is the bound). **Closing step:** after the last gas solve of the time step each record gives the parcel the same correction the gas received, −β_r (φ_new − φ*), so parcels and gas book identical numbers and the pairs cancel to the bit even with one pass. Slopes are capped so that the closing step cannot ask a parcel for more mass or heat than it has; nothing is clipped after the solve. The latent heat is h_v(T) − h_l(T) from one property set (D92, D90 (3)), the evaporated mass carries the vapour enthalpy at the surface temperature, so the energy pair closes too. |
| D251 | **Parcel integrators: exact relaxation and sub-cycling.** Inside a sub-step, with the gas frozen: velocity by the exact solution of the linear drag relaxation with gravity (from SP1), temperature by the exact relaxation towards its equilibrium (the wet-bulb state when evaporating), diameter by the d² law of the sub-step (Godsave 1953; Spalding 1953) with the end of the droplet's life found inside the sub-step, so mass is never negative and the last mass is deposited where it evaporates. Coefficients (Re, Nu, Sh, B_M, properties) are frozen per sub-step and renewed between sub-steps; a sub-step ends at a cell face, at the end of an eddy (D243), or when mass, temperature or Reynolds number have changed by more than a set fraction (case keys, defaults fixed with SP5 from the d² and wet-bulb tests). The flow time step is thus never limited by a response time of a droplet. |
| D252 | **Passes inside the time step.** A parcel pass runs in every outer iteration of the time step (PISO with outer iterations, transient SIMPLEC, coupled solver), always from the committed parcel state of the old time level and with the same random numbers; only the last pass is committed, followed by the closing step. One pass per step is the linearised implicit coupling: stable, bounded, conservative, first order in Δt. Repeated passes remove the linearisation and splitting error: at convergence the step is the simultaneous implicit solution of gas and parcels, whatever the number of passes it took. The *coupling lag* (the largest relative closing correction of the step, per quantity) is the convergence measure, written to the run record (D69). Keys `[spray.coupling] passes` (a number, or `"converged"` with `tolerance` and `max_passes`). **Decided by the maintainer 2026-11-24: the default is `passes = "converged"` with `max_passes = 3`**; the passes of a step end as soon as the coupling lag is below the tolerance, so a weakly coupled step costs one or two. A fixed number (for example 1 for a cheap run) stays a case key; a step that ends at `max_passes` with the lag above the tolerance is counted as a health number. |
| D253 | **Cell limiters as a safety net that is counted.** With D250 and D252 the gas should never pass saturation or the droplet temperature. As a guard for one-pass runs: when the summed evaporation of a cell exceeds what brings the gas of the cell (its mass plus the mass leaving through its faces in the step) to adiabatic saturation, every evaporation record of the cell is scaled by one factor, computed from the sorted sums, applied in the closing step to both books alike. Each action is counted (cells, mass held back) as a health number; a diagnostic code warns when the limiter acts in more than a set share of the spray cells, with the fix "more passes or a smaller step". A limiter never acts silently, and a converged run must show zero actions (V28.5). |
| D254 | **Error checks that ship with the feature.** (1) Closed-box ledger: mass, vapour mass, momentum and enthalpy of gas plus parcels plus walls close to round-off every step, pair residuals zero to the bit (D244). (2) Coupling lag and limiter counts in the run record and the log. (3) Coupling numbers Δt/τ per cell (momentum, heat, vapour) from the summed slopes, their maximum and the share of cells above 1, as health numbers: they say where the step is stiff. (4) A time-step study as a registered tool (D196): the same case at Δt, Δt/2, Δt/4 with the observed order and the estimated error of the monitored quantities (gas temperature at the outlet, evaporated mass, penetration). (5) Liquid volume fraction per cell and parcels per cell as health numbers: the point-source assumption and the statistical error of the sources depend on them (Subramaniam 2013; Garg, Narayanan and Subramaniam 2009). |

### 5.3 The time step with a spray (SP2 onward)

1. Inject the parcels of the step (keyed random numbers, D242, D246).
2. Outer iteration k = 1, 2, …:
   a. freeze the gas state (D248);
   b. parcel pass from the old level for all parcels, in parallel (D247, D251);
   c. sort the records and sum them per cell (D249);
   d. solve momentum, pressure correction (with the mass source), energy,
      vapour and turbulence with the sources in value-plus-slope form (D250);
   e. coupling lag below the tolerance, or the largest number of passes
      reached: leave the loop (D252).
3. Closing step on the records of the last pass; commit the parcels; book the
   ledger pairs; limiter counts and coupling numbers to the run record (D250,
   D253, D254).

Memory: the records of one pass live until the closing step (one record per
parcel sub-step). Cost: one tracking of all parcels per pass; the sort is by
integer keys.

### 5.4 What is promised and what is not

| Property | Promise | How it is shown |
|---|---|---|
| Independent of parcel order and thread count | bitwise, by construction | V28.1 |
| Conservation gas + parcels | round-off, every step, any number of passes | V28.2 |
| Stable and bounded at any Δt | yes (implicit relaxation form) | V28.3, V28.5 |
| Right state for Δt far above the coupling time | the saturated (equilibrium) state, approached without overshoot | V28.5 |
| Order in time of the coupling | first order with the state at the new level (the robust choice: a mid-step state is second order for small Δt/τ but swings around saturation for large Δt/τ); measured, not assumed | V28.4 |
| Second order in time for the coupling | not promised in v1; open point | – |
| Independence of the mesh when droplets are not small against the cells | not promised in v1 (point-source limit); health number, open point | D254 (5) |

### 5.5 Checks (open literature or exact; added to the list of section 3)

| Check | Content | Passes when |
|---|---|---|
| V28.1 | Order independence: the two-way cases with the parcel arrays shuffled, on 1 to 4 threads, CPU and GPU | fields, parcels and ledger bitwise equal on CPU; CPU = GPU as D245 |
| V28.2 | Closed box with a falling, heating droplet cloud: ledger pairs | residuals zero to the bit; accounts to round-off, with 1 and with converged passes |
| V28.3 | Dense cloud at rest in a gas stream, loading Z from 0.1 to 10, Δt/τ_p from 0.1 to 100: common final velocity | momentum balance of the mixture to round-off; no overshoot of the gas velocity past the parcel velocity |
| V28.4 | Time-step study of V28.3 and V28.5 at Δt, Δt/2, Δt/4 | observed order at least 0.9; error at the practical step quoted in the gallery entry |
| V28.5 | Saturation: droplet cloud in a closed adiabatic box with more liquid than the gas can take, Δt/τ_sat from 0.1 to 100 | end state equal to the adiabatic saturation state from the enthalpy balance (exact, own calculation); vapour fraction monotone in time; zero limiter actions with converged passes |
| V28.6 | Two neighbour cells, the cloud in one only: the same result whichever cell is listed first, and on two threads | bitwise |
| V28.7 | d² law and wet-bulb temperature of a single droplet (section 3) at 1 and at 100 sub-steps | diameter history within the band of the exact law for both |
| V28.8 | One pass against converged passes on the spray-cooling case (D240) | the difference of the outlet gas temperature and of the evaporated mass quoted; coupling lag falling with the passes |

### 5.6 Comparison with the open literature

| Approach | Order independent | Stable for Δt above the coupling time | Conserves | Source |
|---|---|---|---|---|
| Parcel loop that updates the cell while it runs | no | partly | yes | the maintainer's experience with earlier codes; not traced to a publication here |
| Frozen gas, lagged explicit sources (PSI-cell), steady or transient | yes | no | yes | Crowe, Sharma and Stock 1977 |
| Implicit coupling of droplet and gas velocity | yes | for momentum; heat and mass as each source treats them (to be read before SP5) | yes | Dukowicz 1980; Amsden, O'Rourke and Butler 1989 |
| Linearised (value plus slope) sources in the finite-volume equations | yes | yes for the linear part | only with a matching parcel update | Patankar 1980 |
| **This design:** frozen state, sorted deposit, relaxation form on all three exchanges, closing step, passes to convergence | yes, bitwise | yes | round-off | D247 to D254 |

State of the art beyond this, noted and not in v1: (a) estimating the sources
with grid-free kernels and controlling the number of parcels per cell, so that
the result converges with the mesh and the parcel count (Garg, Narayanan,
Lakehal and Subramaniam 2007; Garg, Narayanan and Subramaniam 2009;
Subramaniam 2013; diffusion smoothing of Sun and Xiao 2015 is already in
D117); volume-filtered formulations for dense flows (Capecelatro and
Desjardins 2013); (b) the correction for the disturbance a droplet makes in
its own cell, needed when droplets are not small against the cells, as in
large-eddy simulation on fine meshes (Gualtieri et al. 2015; Horwitz and Mani
2016; Ireland and Desjardins 2017; Balachandar, Liu and Lakhote 2019), with
the same question open for the vapour and the temperature a droplet makes
around itself; (c) exponential time integrators for stiff relaxation terms
(Hochbruck and Ostermann 2010), of which D251 uses the simplest member. We
know of no open publication that promises bitwise independence of parcel
order, thread count and partition for a two-way coupled evaporating spray;
that part is our own requirement (D93).

Clean room (D201): everything above is from papers and reports in the open
literature or is our own derivation; no behaviour or constant of a commercial
code is used. Of the citations, Crowe, Sharma and Stock (1977), Garg et al.
(2009), Subramaniam (2013) and Horwitz and Mani (2016) were checked against
their records on 2026-11-24; the others are cited from the literature record.
Each source is checked and read in full before its model is built, and a
number taken from one gets its `DATA_SOURCES.md` entry first.

### 5.7 Open points of this section

- Default number of passes (D252): answered, converged passes with at most 3
  (2026-11-24).
- Steady sprays: answered, own stage SP2s after SP2 (2026-11-24, D260 in
  section 5.9).
- **Second order in time for the coupling:** an exponentially fitted weight
  between the old and the new gas state per cell (exact for the linear
  relaxation at any Δt/τ) is the candidate; decided from V28.4.
- **Mesh independence of the sources** (kernels, parcel number control) and
  the **self-disturbance correction**: with the first case in which droplets
  are not small against the cells.

### 5.9 Steady spray mode (D260, stage SP2s)

Decided by the maintainer 2026-11-24: the steady mode is a stage of its own,
built after the transient two-way coupling of SP2 is verified. The lead
application (D240) has a steady state; marching it in time costs many more
parcel passes than the steady iteration.

| ID | Decision |
|---|---|
| D260 | **Steady spray mode: source-in-cell iteration with the rules of this section.** The method is the steady form of Crowe, Sharma and Stock (1977). (1) A parcel is a *stream*: it stands for a number of droplets per second, is started at its injector and followed through the frozen steady gas field until it leaves, sticks or has evaporated. The records hold the change of the mass flow, the momentum flow and the enthalpy flow of the stream per sub-step; the pass is the same pure function, the deposit the same sorted sum (D247 to D249), the integrators those of D251 with sub-steps along the path instead of a flow time step. (2) A stream pass runs every N outer iterations of the steady solver (key, default fixed with the stage). The gas takes the sources in value-plus-slope form (D250), relaxed between passes, S ← S + α (S_new − S), with α a case key. (3) Random numbers of dispersion are keyed by the stream and the event only, never by the pass: every pass draws the same numbers, so the iteration has a fixed point and can converge. The statistical error of that one sample is estimated by repeating the run with other stream keys (an option of the study tools, D196), not by drawing new numbers each pass. (4) There is no closing step, because no parcel state is carried from pass to pass. The ledger pair is exact for the sources of the last pass; the difference between the relaxed sources the gas holds and those of the last pass is the *source residual*, written to the run record, and part of the convergence test next to the monitored quantity (outlet gas temperature, evaporated mass flow). A run with the source residual above its tolerance is not reported as converged. (5) The limiter of D253 and the health numbers of D254 apply with the residence time of the gas in the cell in place of the time step. (6) Refused with a diagnostic code: injectors with a mass flow that changes in time, and a steady spray together with `[time]`; the wall film (SP7) stays transient. |
| D320 | **Liquid property library.** `[spray.liquid] library` names a liquid whose vapour pressure (Wagner form), density, specific heat, latent heat (Watson form), surface tension and viscosity are functions of T in `foehncraft_kernels::liquids`: water, n-heptane, n-dodecane fitted to the NIST WebBook saturation tables (SRD 69), iso-octane to published correlations; only our coefficients are stored, the generators are in `tools/data`. The vapour species come with it in the species data (NASA-7 plus Lennard-Jones data from the critical point). Density and specific heat stay constant in a run at `property_temperature`. *Implemented 2026-11-26 (PR #321).* |
| D321 | **Film state of the evaporation: the 1/3 rule by default.** The gas properties of the droplet evaporation are those of the film at T_m = T_p + (T_g − T_p)/3, Y_m = Y_s + (Y_g − Y_s)/3 (Hubbard, Denny, Mills 1975; Yuen, Chen 1976), from the species data: Wilke viscosity, Mathur conductivity, mass-weighted c_p of vapour and carrier, ideal-gas density at p₀, Chapman-Enskog diffusivity of the vapour in the carrier with Blanc's law; Re with the film viscosity. The carrier is the start composition without the vapour. `[spray.evaporation] film = "far-field"` keeps the gas properties of the cell. Maintainer's choice of 2026-11-25 (card), after the single-droplet comparison with Chauveau et al. (2019): far field +93 to +99 %, 1/3 rule +39 to +59 %. The high-pressure film of D313 builds on it. *Implemented 2026-11-26 (PR #321).* |

Checks of the steady mode:

| Check | Content | Passes when |
|---|---|---|
| V28.9 | Spray-cooling duct: steady mode against the transient run marched to its steady state | outlet gas temperature and evaporated mass flow agree within the error band of the time-step study (V28.4) |
| V28.10 | Order independence of the steady mode: streams shuffled, 1 to 4 threads | fields and ledger bitwise equal |
| V28.11 | Duct with more water than the gas can take | outlet at the adiabatic saturation state (exact, from the enthalpy balance); vapour fraction monotone along the duct; source residual falling with the passes |

### 5.8 References of this section

- Amsden, A. A., O'Rourke, P. J., Butler, T. D. (1989), KIVA-II: a computer
  program for chemically reactive flows with sprays, Los Alamos report
  LA-11560-MS.
- Balachandar, S., Liu, K., Lakhote, M. (2019), Self-induced velocity
  correction for improved drag estimation in Euler-Lagrange point-particle
  simulations, J. Comput. Phys. 376.
- Capecelatro, J., Desjardins, O. (2013), An Euler-Lagrange strategy for
  simulating particle-laden flows, J. Comput. Phys. 238.
- Crowe, C. T., Sharma, M. P., Stock, D. E. (1977), The particle-source-in
  cell (PSI-CELL) model for gas-droplet flows, J. Fluids Eng. 99, 325–332.
- Crowe, C. T., Schwarzkopf, J. D., Sommerfeld, M., Tsuji, Y. (2011),
  Multiphase Flows with Droplets and Particles, 2nd ed., CRC Press.
- Dukowicz, J. K. (1980), A particle-fluid numerical model for liquid sprays,
  J. Comput. Phys. 35, 229–253.
- Garg, R., Narayanan, C., Lakehal, D., Subramaniam, S. (2007), Accurate
  numerical estimation of interphase momentum transfer in Lagrangian-Eulerian
  simulations of dispersed two-phase flows, Int. J. Multiphase Flow 33.
- Garg, R., Narayanan, C., Subramaniam, S. (2009), A numerically convergent
  Lagrangian-Eulerian simulation method for dispersed two-phase flows, Int. J.
  Multiphase Flow 35.
- Godsave (1953) and Spalding (1953): as in the
  [assessment](../assessments/spray-module-prerequisites.md).
- Gualtieri, P., Picano, F., Sardina, G., Casciola, C. M. (2015), Exact
  regularized point particle method for multiphase flows in the two-way
  coupling regime, J. Fluid Mech. 773.
- Hochbruck, M., Ostermann, A. (2010), Exponential integrators, Acta Numerica
  19.
- Horwitz, J. A. K., Mani, A. (2016), Accurate calculation of Stokes drag for
  point-particle tracking in two-way coupled flows, J. Comput. Phys. 318.
- Ireland, P. J., Desjardins, O. (2017), Improving particle drag predictions
  in Euler-Lagrange simulations with two-way coupling, J. Comput. Phys. 338.
- Patankar, S. V. (1980), Numerical Heat Transfer and Fluid Flow, Hemisphere.
- Subramaniam, S. (2013), Lagrangian-Eulerian methods for multiphase flows,
  Prog. Energy Combust. Sci. 39, 215–245.
- Sun, R., Xiao, H. (2015): as in the assessment.

## 6. Status

SP1 is implemented (PR #242): `[spray]` with an included injector file,
constant droplet size, exact exponential drag and gravity, tracking on the
locator, stick, rebound, escape and mirror boundaries, mass book (SPR001 to
SPR005), parcel output and restart, the parcel kernel in `foehncraft-kernels`
with a PTX entry (the PTX job of CI checks it; it cannot be built without
the nightly toolchain). SP1 is one way: parcel passes as pure functions
(D247) and the sub-step rules of D251 are the first change of SP2.

SP1b is implemented (PR #265): `[[datum]]` with the types plane (cut or
`snap = "faces"`, disc and rectangle bounds), surface (STL), boundary,
interface, box, cylinder and zone; crossing records sorted and summed per
datum (D263), the same bits on 1 to 4 threads and after a restart;
gas reports (mass, volume and momentum flow, area and mass-flow averages,
`cut_error`); moment sums, D10 to D43, Dv10/50/90, histograms and a
Rosin-Rammler fit (D264); penetration, centre of mass and D32 in the flow per
step (D265); `datums.csv`, `spray_hist_<datum>.csv`, `result.json`, the tools
`post.datums`, `post.droplets`, `post.penetration`, the Python reader
and `series.json` (D267, D268). The Rosin-Rammler and log-normal size
distributions of D242 came with SP1b instead of SP2, as the statistics need
them. Not yet: the enthalpy flow (with the energy equation and SP5), size
classes, averaged fields (D266, SP2), `datums.csv` in steady runs, V28.12
(terminal velocity against Gunn and Kinzer).

**SP2a (two-way momentum) is implemented** (`[spray.coupling] mode =
"two-way"`). What the code does, against section 5:

- D247: `Spray::pass` is a pure function of the committed parcels, the
  frozen gas and the keyed random numbers; it can be repeated, and the
  release of the step is booked only at the commit. In two-way mode a
  sub-step ends at every interior face (the walk stops at the first mesh
  face; the time to the face is found by bisection of the displacement
  projected on the chord) and when the particle Reynolds number would move
  by more than `max_reynolds_change` (D251). A parcel that moves along a
  face while the gas of the two cells pushes it back across from both
  sides (the inflow towards the axis of a spray whose nozzle sits on a
  face) would end every sub-step at that face after no time: a sub-step
  that starts on a face and goes straight back through it takes its full
  length, booked to the cell it starts in (found in the gallery run, test
  `a_parcel_on_a_face_pushed_from_both_sides_is_not_lost`).
- D248: the first pass sees the state the step starts from (not
  extrapolated); every later pass the current iterate.
- D249: one record per sub-step (cell, parcel, sub-step, impulse to the gas,
  slope), sorted by (cell, parcel, sub-step) and summed per cell; the
  impulse is the parcel's own change of momentum less its weight.
- D250: S = E − β (u − u*) with β on the momentum diagonal and in D, D^s,
  D^C of the pressure correction, as the porous sink. **Implementation
  note (own derivation):** the slope of record r is the parcel mass times
  (1 − e_r) Π_{l>r} e_l, e = exp(−δt/τ), so the slopes of one parcel add up
  to m (1 − Π e) ≤ m. Summing m (1 − e_r) over the sub-steps, as the plain
  reading of D250 gives, counts the same relaxation several times and lets
  the closing step push a parcel past the gas when it takes many sub-steps
  in one cell. The closing step gives every parcel Σ_r slope_r (u_new − u*)
  of momentum, the gas books the same per-cell numbers with the opposite
  sign, and the ledger pair is formed from one number (zero to the bit).
- D252: passes run inside the outer iterations of PISO (they add outer
  iterations when `solver.outer_iterations` is smaller than the passes
  asked for); default `"converged"` with `max_passes = 3` and `tolerance =
  1e-3`. **The coupling lag** is the largest change of the gas velocity of a
  spray cell against the frozen state of the pass, over the larger of the
  largest slip velocity of the spray cells and the largest change of the
  gas velocity of a spray cell in the step. Both are velocity differences,
  so the measure does not depend on the frame; the slip alone made the lag
  blow up once the parcels had come to rest in the gas (every step then
  ran to `max_passes`). Steps that end at `max_passes` above the
  tolerance are counted (SPR006 warns once).
- D254: (1) the ledger books the gas side on the momentum accounts of the
  domain (which hold only this term for now, so their imbalance is the
  force of the boundaries on the gas), the power of the source on the
  kinetic energy, and accounts `parcel momentum x/y/z` (weight, injected,
  left the flow, walls, from the gas) that close to round-off; (2) passes
  and lag in `spray.csv` and `result.json`; (3) the largest Δt/τ of the gas
  momentum of a spray cell, β Δt / (ρ V), and the share above 1; (5) the
  largest liquid volume fraction and parcels per cell. The time-step study
  as a registered tool (4) is not built yet.
- Checks: V28.1 (shuffled storage, 1 and 4 threads, hex, tet, poly), V28.2
  and V28.3 (a cloud and its gas in one cell, loading 0.1 to 10, Δt/τ_p 0.1
  to 100: one common velocity, the mixture momentum kept, no overshoot, one
  pass and converged passes), V28.4 (observed order 0.99 and 0.96 with Δt =
  2, 1, 0.5 ms against a fourth-order reference). V28.5 to V28.8 belong to
  the heat and vapour exchange (SP2b, SP5).
- Not yet with two-way coupling (SPR003): solid regions, the frozen-flow
  march, algorithms other than PISO.

**SP2b (two-way sensible heat) is implemented** (with `[energy]` in the
case). What the code does, against section 5:

- Parcel temperature (`[liquid] specific_heat`, `temperature`,
  `[[injector]] temperature`): in each sub-step T_p relaxes exactly towards
  the gas temperature of the cell, T_p ← T_g + (T_p − T_g) e^(−δt/τ_h),
  τ_h = ρ_p c_p d² / (6 Nu k_g), Nu = 2 + 0.6 Re^½ Pr^⅓ (Ranz and Marshall
  1952), with the gas conductivity and specific heat of the fluid. One-way
  sprays in a case with `[energy]` carry the temperature too.
- D249, D250: the record of a sub-step holds the heat the gas gives,
  m c_p (T_p − T_p,new), and its slope m c_p (1 − e_r) Π_{l>r} e_l with the
  thermal decays (own derivation, as for momentum: the slopes of a parcel
  add up to m c_p (1 − Π e) ≤ m c_p). The energy equation takes the per-cell
  sums as Q = E − G (T − T*), G on the diagonal (the same path as the
  linearised wall flux), so the gas temperature of a cell stays between its
  old value and the droplets at any step.
- D252: the passes stop when the larger of the velocity lag and the
  temperature lag is below `tolerance`; the temperature lag is the largest
  change of the gas temperature of a spray cell against the frozen state of
  the pass, over the largest temperature difference between gas and
  droplets of the spray cells. **Finding:** the temperature lag falls by a
  factor of about 7 per pass, the velocity lag by about 10; in the
  spray-cooling case most steps end at 1 to 3 × 10⁻³ with the default
  `max_passes = 3`, so the case asks for 4 (one unconverged step in 500).
  The default stays 3 (maintainer decision of 2026-11-24).
- Closing step: each parcel gets Σ_r slope_r (T_new − T*) of heat; the gas
  books the same per-cell numbers. D254 (1): the pair 2804 (`spray heat` on
  the thermal energy account of the energy equation, `from the gas` on the
  account `parcel enthalpy` with injected and left the flow) is zero to the
  bit; two modules that book the same quantity and region share one account,
  so the thermal account of the gas closes with the spray term (2 × 10⁻¹²
  in the gallery run). D254 (3): the largest G Δt / (ρ c_p V) is
  `coupling_heat_max`.
- Checks: V28.1 (shuffled storage and 1 and 4 threads with heat), V28.2
  and V28.3 for heat (a cloud of water at 300 K in one cell of air at 400 K,
  loading 0.1 and 1, Δt/τ_h 0.03 to 23: the mixture temperature to 10⁻⁹, no
  overshoot, enthalpy kept), V28.4 for heat (observed order 1.03 and 1.04
  against a fourth-order reference). V28.16 (one-dimensional dusty-gas
  duct) is not run yet.

**SP2c (eddy-lifetime dispersion) is implemented** (`[spray] dispersion =
"eddy-lifetime"`, D243). What the code does:

- A parcel carries its eddy (fluctuation u′, time left, eddy count; in the
  restart file, older files start new eddies). When the time left is used
  up it enters a new eddy: u′ Gaussian with variance 2k/3 per component
  (Box and Muller from the keyed uniforms of the stream `SPRAY_DISPERSION`,
  key (parcel id, eddy count, component), D121), and the interaction time
  min(τ_e, t_cr) of Gosman and Ioannides (1983), τ_e = l_e/(2k/3)^½,
  l_e = C_μ^¾ k^{3/2}/ε, t_cr = −τ_p ln(1 − l_e/(τ_p |u + u′ − v|)).
  Drag, Reynolds number and heat transfer see u + u′; a sub-step ends at
  the end of the eddy (D251).
- The pass stays a pure function of the committed parcels (D247): every
  pass of a step and every thread count draws the same eddies; the eddy
  the last pass ends in is committed with the parcel.
- k from the turbulence model, ε from the k-ε family or β* k ω (SST,
  β* = 0.09); in 2D cases the eddies have no z component. Without a
  turbulence model the key is refused (SPR001). A cell without turbulence
  (k or ε zero) ends the eddy.
- **Found while testing:** an eddy whose rest is below what the step time
  resolves made its sub-steps stall (the parcel was counted as lost after
  `max_substeps`); such an eddy now ends.
- Checks: tracers in uniform turbulence against the exact random walk of
  the model (within 1.3 % for 20 000 tracers), crossing trajectories (heavy
  droplets change eddies more often and spread less), V28.1 with eddies,
  bitwise equal flow runs on 1 and 4 threads. V28.17 (grid turbulence
  against Snyder and Lumley 1971 and Wells and Stock 1983) is not run yet.

**SP2s (steady spray mode) is implemented** (`[spray.steady]`, D260,
PR #292). What the code does:

- A steady case (no `[time]`) with `[spray.steady]` runs the spray as
  streams: `streams` per injector (default 200), each with the mass flow
  of its injector over the number of streams as its mass (kg/s) and the
  droplets per second as its count. Size, direction and eddies come from
  the keyed random numbers with the key seed · streams + k, the same in
  every pass (D260 (3)); `seed` draws another sample.
- A stream pass is the parcel pass of D247 with `max_time` in place of the
  time step: the stream is followed to the end of its path (escape,
  stick, or `max_time`, SPR007 for unfinished streams). The records hold
  mass, momentum and enthalpy flows; the slopes are the local ones of each
  cell visit, m (1 − Π e), since there is no closing step (D260 (4)).
- The pass runs every `every` outer iterations of SIMPLE or SIMPLEC
  (before the iteration, so the gas solves with fresh sources). The gas
  holds S(u) = A − β u and Q(T) = B − G T, relaxed by
  A ← A + α (A_new − A) and the same for β, B, G; the first pass after a
  start or a restart takes α = 1.
- Source residual: Σ|S_held(u*) − S_new(u*)| / Σ|S_new(u*)| over the cells
  of either, for momentum and heat; the larger is the convergence number.
  A steady run stops only when the gas has converged and the source
  residual is below `tolerance`; when the gas is done first the next pass
  is forced. SPR008 marks a run that ends above it as not converged.
- After the loop a final pass gives the books, the ledger and
  `spray/streams.vtkhdf`. The ledger books the pairs 2801 to 2803 and 2804
  with the sources of the final pass on both sides (zero to the bit) and
  the difference to the held sources as the term `spray source residual`
  on the gas accounts; the parcel accounts are in the region `parcels`.
  The rows go into the rates before the steady balance, so the kinetic
  energy view sees the spray power.
- D254 health numbers in the steady form: β/ṁ and G/(ṁ c_p) with ṁ the gas
  mass flow through the cell (half the sum of the face mass flows), and
  the liquid fraction from the volume flow times the dwell time.
- Checks: V28.10 (spray.csv and residuals bitwise equal on 1 and 4
  threads, with and without dispersion), a single stream in a uniform gas
  against a fourth-order integration (residence time to 7 × 10⁻⁵), another
  seed draws other streams, the refusals of D260 (6). V28.9 without
  evaporation: the steady duct takes 6.39 kW from the air against a time
  mean of 6.58 ± 0.03 kW of the transient run after 0.3 s (2.9 % apart,
  more than the error band of the transient run; the transient flow does
  not come to rest behind the nozzle). Not met in the strict form; a gap of
  the feature row. V28.11 waits for evaporation (SP5).
- Not yet: evaporation in the steady mode (SP5), liquid at datums, solid
  regions, the coupled algorithm, the seed study as a tool (D196), several
  processes.

**SP4 and SP5 (evaporation) are implemented** (`[spray.evaporation]`,
`[liquid] vapour`, PR #306). What the code does, against section 5:

- SP4: the liquid names the species it evaporates into (`vapour`, which
  must be a species of `[species]`; refused otherwise). Water only: the
  saturation pressure is the IAPWS-95 auxiliary equation (Wagner and Pruss
  1993), the latent heat and the vapour enthalpy follow from it and the
  specific heats. The gas density stays the low-Mach density of the case;
  the examples use constant density (see the gaps).
- Models: `ranz-marshall` (Sh = 2 + 0.6 Re^½ Sc^⅓, the Spalding number
  B_M) and `abramzon-sirignano` (film corrections of the Nusselt and
  Sherwood numbers, B_T from B_M). In each sub-step the diameter and the
  temperature advance with frozen coefficients; the temperature relaxes
  exactly towards the wet-bulb balance of the sub-step and never passes it.
  D251: a sub-step also ends when d² would fall by more than
  `max_mass_change` or T_p change by more than `max_temperature_change`; a
  droplet below `min_diameter` whose life ends in the sub-step is gone.
- D249, D250: the record of a sub-step carries the vapour mass, its slope
  (the decay share of the vapour against Y − Y*), the momentum and the heat;
  the vapour enters the continuity equation, the species `vapour` and the
  energy equation (with the vapour enthalpy) as S = E − β (Y − Y*), β on
  the diagonal of the species row. A negative explicit part (condensation)
  goes on the diagonal too (Patankar), so Y stays non-negative. **Finding
  (own derivation):** the enthalpy that goes with the slope and with the
  closing step must be the enthalpy the droplet loses with that vapour,
  h_l + e_h (h_v − h_l) with e_h the thermal decay of the record (c_l T_p in
  the steady mode). With the vapour enthalpy h_v the gas of a saturating
  cell ended about 2 K below the adiabatic saturation; with the pure liquid
  enthalpy the droplet temperature oscillated by ±0.7 K from step to step.
- D253 limiter: the vapour a cell can take in a step is its own capacity to
  the droplet saturation (ρ V (Y_sat − Y)) plus that of the gas flowing in
  over the step (inflow mass times its capacity); in the steady mode only
  the inflow per second. Where the records ask for more, every record of
  the cell is scaled by one factor; the held-back mass stays liquid at the
  enthalpy of the record, so the books close. The transient mode revives
  droplets that would have vanished; the steady mode books the held-back
  flow as `held back (D253)` on the parcel mass account. SPR009 warns.
- D252: the coupling lag includes the vapour lag. **Finding:** with
  evaporation the lag falls by a factor of only about 1.6 per pass, so the
  evaporating duct asks for `max_passes = 8` (24 of 200 steps still end
  above 10⁻³).
- D254 (1): the vapour pair 2805 (`spray vapour` on the gas mass account,
  `to the gas` on `parcel mass`) is zero to the bit, as 2801 to 2804.
- Steady mode (SP2s): the stream pass carries the evaporation; vapour and
  its enthalpy are held in value-plus-slope form and relaxed as momentum and
  heat.
- Checks:
  - V28.1 with evaporation (spray.csv bitwise on 1 and 4 threads) and a
    crowded cell: 4000 parcels evaporating in one cell at 95 % saturation
    never push the vapour past the droplet saturation nor below zero, and
    sources and parcels are bitwise equal on 1 and 4 threads with the
    storage reversed (the records are sorted by cell and parcel before they
    are summed, D249; each parcel sees only the frozen gas, D247).
  - V28.2/V28.5 in the form of a closed box: water and energy kept to
    10⁻¹², the gas at the exact adiabatic saturation within 0.003 K.
  - V28.7: a 100 µm droplet in dry air at 500 K follows the d² law (local K
    within 1 %) at the wet-bulb temperature.
  - V28.11 (steady saturation duct, case `spray_saturation_steady`)
    **partly met**: the air saturates within about 5 cm and the vapour
    stays flat to the outlet (no overshoot); the outlet is 309.3 K and
    Y = 0.0375 against the exact adiabatic saturation 310.1 K and 0.0395
    (own ODE c_p dT = −L(T) dY/(1 − Y)). The source residual falls from 1
    to about 0.1 and stalls, also with relaxation 0.1 and 0.03 (SPR008),
    and the limiter still acts in the final pass (SPR009).
  - V28.3, V28.4, V28.6 and V28.8 are not run as named checks with
    evaporation.
- Not yet: other liquids, boiling, the GPU entry of the evaporation kernel,
  several processes; the heat balance closes to round-off only with
  `[species] properties = "fluid"`.

## 7. Breakup (SP6, D300 to D309)

The owner asked on 2026-11-25 for the current best breakup model, for the
ECN Spray G validation (iso-octane, condition G1) and then Spray A. The
assessment of the 2021 to 2025 literature is
`foehncraft/assessments/spray-breakup-assessment.md` (project files): the
Spray G and Spray A studies of Argonne, Sandia and others still use KH-RT
with blob injection; what moved forward is the injector input (measured rate,
area contraction, plume direction and cone) and the distortion of the drops.

| No. | Decision |
|---|---|
| D300 | **Model family: KH-RT with blob injection and a breakup length, plus TAB distortion with dynamic drag** (Reitz and Diwakar 1987; Reitz 1987; Patterson and Reitz 1998; Beale and Reitz 1999; O'Rourke and Amsden 1987; Liu, Mather and Reitz 1993). Keys `[spray.breakup] model = "kh-rt"`, `drag = "dynamic"`; the constants are case keys with the defaults of the model papers (B₀ = 0.61, B₁ = 10, C_RT = 1, C_τ = 1, C_bl = B₁/2), and a validation case calibrates B₁ against its own measurements, not against other codes' tuned values (D201). The surface tension and viscosity of the liquid are case keys (constant) until the liquid property library takes over. Put to the owner on a card on 2026-11-25 with this option recommended; built on it while the answer was open. **The owner chose KH-ACT** (22:13 UTC): D307 adds the nozzle turbulence and the cavitation on top of this core, which stays as described here and as the model `kh-rt`. |
| D301 | **Kelvin-Helmholtz stripping, exact per sub-step.** The wave of Reitz's (1987) curve fits (Λ, Ω from We_g, the Ohnesorge number and the Taylor number) sets the child radius r_c = B₀Λ and the time τ = 3.726 B₁ a/(ΛΩ); with the gas frozen over a sub-step, da/dt = −(a − r_c)/τ is solved exactly, a = r_c + (a − r_c) e^(−Δt/τ). Own check of the fits against the inviscid dispersion relation (Bessel functions): the wavelength within 5 % for every We_g; the growth rate within 3 % below We_g = 0.1 and above 50, up to 46 % high around We_g = 2 (the fit is the sum of the two limits). |
| D302 | **Child parcels from blobs, identities from the parent, mass exact at every split.** A blob keeps its liquid while its KH radius shrinks (the ratio `kh` of the parcel state); when the stripped liquid reaches `shed_fraction` (default 0.03, Patterson and Reitz 1998) of the blob's injected liquid it leaves as a child parcel of drops of radius r_c with the blob's position, velocity, temperature and eddy. A blob strips at most that share in one sub-step (a sub-step limit), so each child carries one to two shares. Children and the drops of other injectors strip in place (smaller drops, more of them). Beyond the breakup length L_b = C_bl (ρ_l/ρ_g)^½ d_hole from the injector the Rayleigh-Taylor wave of the drag deceleration g_t = \|u − v\|/τ_p competes: while C_RT Λ_RT is smaller than the drop it grows, and after C_τ/Ω_RT the drop breaks into drops of the size C_RT Λ_RT (the blob becomes a drop). A child's identity holds the injector (top 16 bits), bit 47, the child index (15 bits) and the parent's injection index (32 bits), so the dispersion and every other keyed random number of a child is its own, and a blob and its children are tracked in the same pure step (D247): the results are the same bits for any thread count and parcel order. |
| D303 | **TAB distortion and dynamic drag.** Every parcel carries y and ẏ of the Taylor analogy (O'Rourke and Amsden 1987: C_F = 1/3, C_k = 8, C_d = 5, C_b = 1/2), solved exactly per sub-step (under-, critically and over-damped); y is held at 1 once reached (no TAB breakup). With `drag = "dynamic"` the drag coefficient is C_D,sphere (1 + 2.632 y), 0 ≤ y ≤ 1 (Liu, Mather and Reitz 1993), in the drag relaxation, the eddy interaction time and the RT deceleration. The Spray A work adds the corrected distortion of Nguyen, Dahms, Pickett and Tagliante (IJMF 2022) on the same y. |
| D304 | **The state in the restart file and the output.** The breakup state (y, ẏ, `kh`, the RT time, the children shed, the blob flag, the injected liquid) is part of the parcel state of the restart file (format `foehncraft.spray-state` 1; a file without it starts undistorted drops); the parcel file has the field `distortion` with the dynamic drag; `spray.csv` counts `child_parcels` and `rt_breakups`. |
| D305 | **Refused in this version:** breakup and dynamic drag in the steady spray mode (a stream cannot shed child streams yet); more than 2³² parcels per injector with breakup (the bits of the identity). |
| D306 | **Injector keys for nozzles.** `hole_diameter` and `area_contraction` C_a: without `diameter` and `size` the parcels start as blobs of √C_a d_hole, without `speed` they leave at U = ṁ/(ρ_l C_a A_hole); `rate` is the mass flow over time, a table [t, ṁ] linear in between, and each parcel carries the exact integral of the table over its share of the injection time; `flow_share` is the share of the injector's liquid in the domain (a quarter or sector model of a jet), with the speed from the whole flow. |
| D307 | **KH-ACT: the nozzle turbulence and the cavitation strip the blobs** (Som and Aggarwal, Combust. Flame 157, 2010, eqs. 7 to 19; formulation also in Magnotti's thesis, Georgia Tech 2017, OSTI 1785739, eqs. 1.13 to 1.25). Key `model = "kh-act"`. A blob is stripped by the mechanism with the largest ratio of length to time scale at the start of each sub-step: Kelvin-Helmholtz (a − r_c)/τ_KH as in D301; the turbulence the liquid brings out of the hole, L_T/τ_T with L_T = C_μ k^1.5/ε, τ_T = C_μ k/ε (C_μ = 0.09); the cavitation R_CAV/τ_CAV with R_CAV = r_hole (1 − C_a)^½, τ_CAV the shorter of the Rayleigh collapse time 0.9145 R_CAV (ρ_l/p_v)^½ and the burst time (r_hole − R_CAV)/(2k/3)^½. Turbulence and cavitation take the radius down at da/dt = −C_T,CAV L/τ (key `c_tcav`, default 1) to no less than their child radius (L_T, R_CAV), and the children of the stripped liquid have that radius (KH children keep B₀Λ). The turbulence of a parcel decays without production, convection or diffusion as in the standard k-ε model (C_ε2 = 1.92); our own solution of those two equations is exact, k = k₀ (1 + t/T)^(−1/(C_ε2−1)), T = k₀/((C_ε2 − 1) ε₀), so L_T grows as (1 + t/T)^0.457, τ_T = τ_T0 + 0.0828 t, L_T/τ_T = √k, and the stripped radius over a sub-step is the exact integral of √k (the papers' constants 0.0828 and 0.457 follow). k₀ and ε₀ at the exit: `nozzle_k`, `nozzle_epsilon` (for example from a nozzle flow simulation, as Som and Aggarwal used), or the force balance of Huh and Gosman (1991, as given by Magnotti eqs. 1.22 and 1.23), k₀ = U²/(8 L/d) B, ε₀ = 0.27 U³/(2L) B, B = 1/C_d² − K_c − (1 − s²), from `hole_length`, `discharge_coefficient` and `inlet_loss` K_c (default 0.5, the textbook loss of a sharp-edged entrance; our choice, the papers give no value), with U the blob's speed from the hole and s → 0 (the hole fed from a much larger sac, our reading). `[spray.liquid] vapour_pressure` sets the collapse; without it the cavities only burst. Refused: C_d above C_a, losses that leave no friction (B ≤ 0), half of a pair of keys. Only blobs feel the turbulence and the cavitation (Som and Aggarwal: primary breakup acts on the parent parcels); children and drops follow KH-RT. `spray.csv` counts `act_children`, the children the turbulence or the cavitation stripped. |
| D308 | **The eddy is mirrored with the parcel.** At a mirror boundary (symmetry plane, slip wall) the normal component of the parcel's eddy fluctuation changes sign with its velocity: the parcel goes on in the mirror image of the flow. Before, the fluctuation drove a small drop back into the plane at once and it bounced in tiny sub-steps for the eddy's life: in the gallery spray, injected on the corner of two symmetry planes, 47 million rebounds in 0.3 ms; with the mirrored eddy 25 thousand, the run 8 times faster and the penetration the same to 0.1 mm. |
| D309 | **ECN Spray G, condition G1, on a sector of the injector.** A sector mesh: `[mesh.generate] wedge` with `dimension = 3` revolves a 3D box of hexahedra about the x axis into a sector of the angle; the z intervals are the angular cells, `zmin` and `zmax` the side planes, ordinary boundaries (refused: tetrahedra, which collapse on the axis, and `perturb`, which would move the nodes on the axis apart). The G1 case computes one eighth of the 8-hole injector: a 22.5° sector with the plume of hole 1 on `zmin` (half a plume, `flow_share = 0.5`) and the plane half way to the next plume as `zmax`, both mirror planes of the 8 plumes. **A thin slip core on the axis** (0.3 mm, `ymin` symmetry) with the cells next to it graded: the prisms that meet on the axis have faces of skewness 1/(4Δφ), 5 with 8 angular cells, which the mesh check remediates to first order; in our trials (8 angular cells) they were the first cells to go unstable once the plumes drew gas to the axis, after 0.13 ms (12 angular cells: sooner), and the parcels in the diverging gas then spent `max_substeps` each and filled the memory with their exchange records; with the core the same runs are clean. By the symmetry of the plumes the gas on the axis moves only along it, which the core imposes. **The ECN metrics** (`tools/plot/spray_g_metrics.py`), both along the injector axis from the tip: the liquid penetration is the farthest distance where the projected liquid volume of the 8 plumes (each parcel mirrored and turned by k·45°), seen perpendicular to the plane of plumes 1 and 5, exceeds 0.2e-3 or 2e-3 mm³/mm² (pixels of 0.1 mm and a Gaussian blur of 0.1 mm, our choices: without the blur the penetration moves by less than 0.1 mm between 0.05 and 0.2 mm pixels); the vapour penetration is the farthest cell whose vapour mass fraction exceeds 0.1 % (the ECN model definition, our threshold). **The injector from the measurements:** the ECN rate of injection (our 17-point fit), the plume direction 33° (our mean of the ECN tomography), the full cone 25°, C_a = 0.65 (measured range 0.635 to 0.687); C_d = 0.51 is ours from the measured rate. KH-ACT with the default constants (B₁ = 10), not tuned to the measurements. |

### 7.1 Checks

- Kernel (`foehncraft-kernels`, module `breakup`): the KH fits against the
  dispersion relation (D301); Rayleigh's capillary limit at U = 0; the KH
  radius relaxation exact (two half steps = one step); the RT wave against
  its formulae; the TAB solution against fourth-order Runge-Kutta in the
  under- and over-damped regimes and at critical damping; the steady
  distortion We_g/12; the dynamic drag factor bounded.
- Parcels (`foehncraft-particles/tests/breakup.rs`): blobs shed children of
  the KH size, each of one to three shares, and the liquid in the flow plus
  out is the injected liquid to 10⁻¹³; bitwise equal on 1 and 4 threads with
  the storage reversed, with children, RT breakups and dynamic drag; a
  restart continues with the same bits; the KH radius of a blob follows the
  breakup equation integrated independently by RK4 within 2 % of the
  stripped part over more than half of its life; RT only beyond the breakup
  length; dynamic drag slows the spray; a rate table gives each parcel its
  exact share and the speed ṁ/(ρ_l C_a A).
- KH-ACT, kernel: the decay of the nozzle turbulence against RK4 of the k-ε
  equations without production (k, ε, L_T, τ_T to 10⁻⁹), the integral of
  √k against Simpson's rule, L_T/τ_T = √k; the Huh-Gosman estimate by hand;
  the cavity radius, collapse and burst times; the strongest mechanism
  wins, the radius stops at the child radius, C_T,CAV scales the rate.
- KH-ACT, parcels: a blob stripped by the turbulence alone has the radius
  a₀ − C_T,CAV ∫√k dt to 10⁻¹² of a₀ whatever the sub-steps, its children
  are at least L_T0 and smaller than the blob; with the Huh-Gosman
  turbulence and cavities the liquid is kept to 10⁻¹³, the results are the
  same bits on 1 and 4 threads with the storage reversed, and the blobs
  hold less liquid than with KH-RT; the nozzle keys are checked.
- The gallery case `spray_breakup` (KH-ACT, k-ε) against the penetration
  correlation of Hiroyasu and Arai (1990): 97 % of the liquid at 50.6 mm
  at 1 ms (correlation 51.1 mm), 47.0 mm with KH-RT; up to 24 % ahead at
  0.2 to 0.4 ms. With the SST model of the first runs the turbulence of the
  thin gas jet did not build up in the still chamber (k near its start
  value) and the tip ran 80 % ahead; the case uses k-ε. ECN Spray G (G1)
  is the next part (D309).

### 7.2 Not yet

Collision and coalescence (the drops of the gallery case have a Sauter
mean diameter near 2 µm, small next to measured diesel sprays); a size distribution of the children (the Jackiw and Ashgriz
fragment model for slow drops is a candidate); breakup in the steady mode;
temperature-dependent surface tension and viscosity (with the liquid
property library); radial velocity of the children (they leave with the
blob's velocity); GPU entry of the breakup kernel.
