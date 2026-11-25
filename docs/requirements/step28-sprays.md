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
| SP2s | Steady spray mode (D260): parcel streams through the steady gas field, relaxed sources, checks V28.9 to V28.11. Momentum and sensible heat first; evaporation in the steady mode when SP5 is merged | SP2 merged and V28.1 to V28.4 passed |
| SP6 | Breakup and collision | decision on regime |
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
- SP6 (which atomisation regime first) is undecided.
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
