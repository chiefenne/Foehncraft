# Step 28: Lagrangian sprays, version 1 (L12), specification SP0

Decisions of the maintainer, 2026-11-24, from the assessment
[spray-module-prerequisites](../assessments/spray-module-prerequisites.md)
(PR #228). Specification only; coding of each part starts on the maintainer's
explicit go. Clean-room rules (D201) apply: models from the open literature,
constants as case keys with the literature value as default, no numbers copied
from commercial documentation. Reacting flows are not part of this step.

## 1. Decisions

| ID | Decision |
|---|---|
| D240 | **Lead application: water-spray cooling of a hot gas stream.** It decides the validation list and the first gallery case: injectors spraying water into a hot air duct, gas cooled by sensible heat first (SP2), then by evaporation (SP5). Urea-water for SCR, fuel sprays, painting and scrubbers follow with their own cases. |
| D241 | **Serial release first.** SP0 to SP2, SP4 and SP5 run on one process (threads allowed). A spray case started on more than one rank is refused with a diagnostic code. Parcel migration between ranks (SP3) is planned after MPI W3 to W5 are merged. All parcel kernels are pure functions of (parcel, gas state), so nothing is rewritten for SP3. |
| D242 | **Droplet size input.** Three forms as case keys: constant diameter, Rosin-Rammler (Rosin and Rammler 1933) and log-normal; a table form is not in v1. SP1 ships constant diameter, SP2 adds the two distributions. Sizes are drawn from the counter-based generator (stream `SPRAY_INJECTION`, D121). |
| D243 | **Turbulent dispersion: eddy-lifetime model first** (Gosman and Ioannides 1983), using k and epsilon (or omega) at the parcel, stream `SPRAY_DISPERSION`. **Answered by the maintainer 2026-11-24: both models as values of one case key** (`dispersion = "eddy-lifetime"` or `"random-walk"`, same stream): eddy lifetime ships in SP2, the continuous random walk follows behind the same key. |
| D244 | **The spray ledger pair waits for the ledger as data (D119, stage 7 part 7.0a).** SP2 (two-way coupling) starts only when 7.0a is on `main`; the gas-parcel exchange (momentum, mass, enthalpy) is then booked as ledger rows with opposite signs on both sides. Consequence: 7.0a becomes a prerequisite of SP2, not of SP1. |
| D245 | **GPU kernels from SP1.** Parcel kernels are single-source from the first pull request, compile to PTX in the existing CI job, and have a CPU = GPU check (to the tolerance of the exact-sum policy of 2026-11-23). SP1 is therefore larger than planned; the separate stage SP9 reduces to multi-GPU and performance work. |
| D246 | **Injector keys in an included file** (D114). The case file names the spray file (`[spray] include = "spray.toml"`); injectors (`[[injector]]`: position, direction or cone, mass flow over time, size form, liquid, temperature) and size or rate tables live there. The resolved case (`case.resolved.toml`) holds the merged keys, so a run stays reproducible and agent-readable. |

## 2. Stages after the decisions

| Stage | Content | Starts when |
|---|---|---|
| SP1 | Particle layer, one rank, one-way: SoA parcels, `[[injector]]`, tracker on `Locator`, exponential sub-stepping, Schiller-Naumann drag, gravity, stick and rebound walls, parcel block in the VTKHDF output and restart, **GPU kernels with CPU = GPU check**. Gallery example, case, recipe. | go |
| 7.0a | Ledger as data (D119), in stage 7 | before SP2 |
| SP2 | Two-way momentum and sensible heat, size distributions, eddy-lifetime dispersion, ledger rows for the gas-parcel pair | SP1 and 7.0a |
| SP4 | Minimal vapour species on the low-Mach density, liquid property set for water first | SP2 |
| SP5 | Evaporation (Ranz-Marshall, Abramzon-Sirignano), closed-box mass and energy ledger test, the spray-cooling case | SP4 |
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
it appears in a table.

## 4. Open points

- D243 (dispersion model) is answered: both models as one key (2026-11-24).
- SP6 (which atomisation regime first) is undecided.
- Liquids beyond water (ethanol, n-heptane, n-dodecane) need source entries
  in `DATA_SOURCES.md` before use.
