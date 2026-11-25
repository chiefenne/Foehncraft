# Step 28: spray validation, statistics and viewing (D261 to D268)

Decisions of the maintainer, 2026-11-24, on three questions: are the planned
checks enough for the first spray stages, what statistics must the solver write
(for example the liquid mass flow through a plane or any datum), and what must
MurmelCraft show. Companion of the [spray specification](step28-sprays.md)
(D240 to D260). Specification only; coding of each part starts on the
maintainer's explicit go. Clean-room rules (D201) apply: open literature,
exact solutions and our own derivations only; no number of a commercial code.

## 1. Where we stand (2026-11-24)

- **Checks.** SP1 is covered (Stokes and Schiller-Naumann terminal velocity,
  exponential relaxation, tracker against reference trajectories on four mesh
  types, bitwise on 1 to 4 threads, restart). For SP2 onward the specification
  has V28.1 to V28.11. They are exact or energy-balance checks. Three gaps:
  (a) no measured data set is chosen, and the reference data sets of section 3
  of the specification have no `DATA_SOURCES.md` entry; (b) the lead
  application (D240, water-spray cooling of a hot gas stream) has only the
  energy balance of the duct; we found no open measurements of a water spray
  into hot air; (c) no check compares a two-way coupled flow with
  measurements.
- **Statistics.** The flow solver writes probes, forces, wall quantities,
  balances and running means (D43, D71, D89), but **no integral over a plane,
  surface or region** of the flow: no mass flow, enthalpy flow or area average
  through a section that is not a boundary. SP1 writes the global mass book
  (`spray.csv`) only.
- **Viewer.** MurmelCraft opens a mesh and the fields of one result file
  (slices, iso-surfaces, arrows, probes). It has no parcel points, no stepping
  through a series of result files, and no charts of CSV files.

## 2. Decisions

| ID | Decision |
|---|---|
| D261 | **Validation matrix and data policy** (section 3). Every stage has exact checks, checks against our own derivations (ODE references) and, where open data exist, one measured anchor. **Decided 2026-11-24:** the measured anchor of the evaporating spray is an open set with a liquid other than water (candidate: Sommerfeld and Qiu 1998, read first), plus exact and own references for water. No water-in-hot-air measurements are known to us; if the maintainer supplies some, they are added as a further check. A reference data set enters a table only after it is read in full and has its `DATA_SOURCES.md` entry (D201); its liquid needs a property set with an entry (specification section 4). |
| D262 | **One generic datum for gas and spray.** A named `[[datum]]` is a geometric object that reports can name: a plane (point and normal, optionally bounded as disc, rectangle or by a region of D73), a boundary, an internal face zone, a surface file (STL, with the ray casting of D108 part 2), or a region (D73) for volume statistics. Reports (`[[report]]`, D263 to D265) name datums and quantities. Two kinds of datum surface: **face zones and boundaries are exact** (the sum of the faces, conservative, bitwise like every other sum); **free planes and surface files are cut through the cells**, values from the cell value and its least-squares gradient (D24), with the cut error quoted in the report: it states the plane's mass flow against the flow through the boundaries and the largest relative difference of a plane that coincides with a face zone (V28.13). A `snap = "faces"` key moves a free plane to the nearest face layer of the mesh when the maintainer wants the exact kind. Gas reports: mass flow, vapour mass flow (convective and diffusive part), enthalpy flow, momentum flow, area average, mass-flow-weighted average, and with D90 the entropy flow. |
| D263 | **Liquid through a datum: crossings, cumulative first.** Every sub-step that crosses a datum writes one crossing record (datum, parcel global id, sub-step number, direction, time, position, velocity, diameter, droplets represented, mass, temperature, injector). Position and time of a crossing come from the closed-form path of the sub-step (D251) where the kernel has it, else from the chord (shorter than a cell: sub-steps end at cell faces, D247). Records are sorted by (datum, parcel id, sub-step) and summed in that order (D249): no atomics, bitwise on any thread count and CPU = GPU. **The primary result is the cumulative mass, droplet number and droplet enthalpy through the datum, signed, per direction**; the flow rates are differences over the report interval. The report writes the number of crossings behind every rate, and a relative standard error from batch means (D71): a rate from few crossings says so. Steady mode (D260): a stream carries its own mass flow, the flow through a datum is the sum of the stream flows crossing it and there is no sampling noise from time. Crossings of walls (stick, rebound) are booked by boundary the same way. Mass of liquid, vapour and (with SP7) film through the datums of a duct gives the water balance between two datums (V28.14). |
| D264 | **Droplet statistics: from moment sums.** At a datum (flux-weighted) and in a region (D73, volume-weighted, per time or time-averaged) the report holds the sums Σ N dᵏ for k = 0 to 4 and the mass, velocity and temperature sums, so that the diameters D10, D20, D30, D32 (Sauter), D43, the volume diameters Dv10, Dv50, Dv90 and the span, the mean velocity and its standard deviation per component and the mean temperature follow from sums that are exact reductions (D31), and merge across ranks as sums (SP3). Histograms: number- and mass-weighted, linear or logarithmic bins, bin range and count as keys, with an overflow bin that is counted; cumulative distribution and the Rosin-Rammler fit (D242) of the histogram, with the fit's own residual. Size classes (a list of bounds) give the mass, velocity and temperature of each class separately, as the phase-Doppler measurements of V28.19 are reported. |
| D265 | **Spray-wide numbers per step.** In `spray.csv` and the run record (D69), per step and per injector: parcels and droplets in the flow, injected, evaporated (SP5), stuck, left, lost (the mass book of SP1 stays), evaporation rate (kg/s, from the ledger rows of D244 so it is the same number the energy balance uses), liquid mass in the flow and its centre of mass, **liquid penetration** (the distance from the injector, along the injector axis or radially as a key, within which a set fraction of the liquid mass lies: default 95 %; the farthest parcel is reported too, it is outlier-sensitive and not the primary number), **vapour penetration** (distance to the farthest cell with vapour fraction above a threshold, key, SP4), **evaporation length** (distance within which a set fraction of the injected mass has evaporated, SP5), and the health numbers of D254 (parcels per cell, liquid volume fraction, coupling numbers). Definitions are written in the tool description and in the file header, so that a number can be compared with a published one only when the definition is the same. |
| D266 | **Time-averaged and sampled spray fields.** Parcels are *sampled* into the gas cells at every sub-step, weighted by the time spent in the cell (sorted sums, D249; the sampling is separate from the exchange deposit of the coupling, and costs nothing when no field is asked for). With the averaging window of `[output.statistics]` (D71) the result is written to `statistics.vtkhdf` beside the running means of the gas: liquid volume fraction, droplet number density, mean droplet velocity (mass-weighted), D32, mean droplet temperature, and with SP2 and SP5 the mean momentum, heat and evaporation sources the gas received. The confidence interval of the averaged report quantities follows D71 (batch means); fields carry the number of samples per cell, so a cell with too few parcels is visible (D254 (5)). |
| D267 | **Structured, scriptable output; no interactive-only path.** Written by the run: `spray.csv` (D265), `datums.csv` (one row per report interval, one column per datum and quantity, named `<datum>.<quantity>`), `spray_hist_<datum>.csv` (histograms, one row per bin and interval), the same numbers in `result.json` (section `spray` and `datums`, with the definitions) and as rows of the run record (D69). Tools of the registry (D196): `datum.report` (any datum and quantity from a result or a restart, without a rerun, for gas fields), `spray.statistics` (statistics of a parcel file, from the same sums as the run), `spray.penetration`, and a Python reader in `foehncraft` that returns the tables as arrays. A script or tool that misses a field, datum or quantity lists the ones that exist (maintainer's rule of 2026-11-24). All of it reachable through the MCP server with the roots of D198. |
| D268 | **MurmelCraft: parcel viewer and time player first.** First milestone in the MurmelCraft repository: (1) read the parcel block of the VTKHDF files of SP1 (points with diameter, velocity, temperature, parcel id, injector) and show the parcels as points or glyphs coloured by a parcel field with its own colour bar and size scaling, together with the mesh and the gas fields, as a model-tree item with eye and filters (injector, diameter range); (2) **step through a series of result and parcel files in time**: a player (play, pause, step, loop, speed), a time slider and the time of the shown files in the status bar, with the files of a run found from a manifest `series.json` (file, time, step) that the solver writes beside the output (added with SP1b if it is not there; the viewer also accepts a folder and sorts by the time stored in the files), as tools (`series.open`, `series.go`, `series.play`), usable by the MCP server and in `export.movie`, so that a spray or any transient result can be recorded as a movie. Second milestone, after the maintainer's go: charts and tables of `datums.csv`, `spray.csv` and the histograms (the datum shown in the model tree as a plane with its cumulative mass and rate), and a statistics panel. Contract on our side: file layouts are documented in `docs/` of this repository and versioned (D80); MurmelCraft never needs the solver to be present. |

## 3. Validation matrix (D261)

All references are exact, our own derivation, or open literature. Citations are
from the literature record and **not yet read**; each is read in full and gets
its `DATA_SOURCES.md` entry before any number enters a table or a band.
New checks continue the numbering of the specification (V28.1 to V28.11).

| Check | Stage | Content | Reference | Passes when |
|---|---|---|---|---|
| V28.12 | SP1b | Terminal velocity of water drops in still air, from 0.1 to about 5 mm | measured: Gunn and Kinzer (1949); our own ODE with the correlation of Schiller and Naumann; Stokes limit exact | within the scatter band stated in the paper for the range where the correlation is valid; larger drops (shape) are shown and not claimed |
| V28.13 | SP1b | Datum contract: a plane on the faces of a mesh against the face-zone datum; a plane shifted by one cell; a rotated plane; gas mass flow through each | exact (face zone); closure against the boundary flows | face zone equal to the plane at the faces to round-off; liquid cumulative mass through every plane of a non-evaporating spray equal within one parcel mass; gas mass flow of cut planes within the quoted cut error |
| V28.14 | SP1b, SP5 | Injector statistics: sizes drawn against the input distribution (constant, Rosin-Rammler, log-normal); mass flow through a datum at the injector against the injector's mass flow; with SP5, liquid plus vapour through two datums against the evaporation of the ledger between them | exact (distribution moments, own derivation) | moments within the sampling error; cumulative mass within one parcel mass; water balance to round-off of the ledger |
| V28.15 | SP1b | Statistics from sums: D10 to D43, Dv50, histograms from the parcel file against a direct calculation; 1 to 4 threads, CPU = GPU, after restart | exact | equal to round-off; bitwise across threads |
| V28.16 | SP2 | One-dimensional dusty-gas duct (uniform droplets, constant exchange coefficients, momentum and heat): ODE solution at 1e-10, compared with the CFD duct's mass-flow-weighted averages at datums along the duct | own derivation (Marble-type relaxation, textbook: Crowe, Schwarzkopf, Sommerfeld and Tsuji 2011) | within the error band of the time-step study (V28.4); order independence (V28.1) |
| V28.17 | SP2 | Dispersion in grid turbulence: spread and velocity autocorrelation of particles of several Stokes numbers, with the eddy-lifetime model and the random walk (D243) | measured: Snyder and Lumley (1971), Wells and Stock (1983); the theory of Csanady (1963) for the limits | spread within the scatter of the data for the particle classes whose properties the papers give; the two models differ by a quoted amount |
| V28.18 | SP2 | Two-way coupled particle-laden flow with measurements: a particle-laden round jet or a particle-laden channel; mean and fluctuating velocity of both phases | candidates, to be read and chosen with the maintainer when SP2 starts: Lau and Nathan (2014) for a jet; Kulick, Fessler and Eaton (1994) for a channel. **Our model has no turbulence modulation by particles** (open point, section 5): comparisons are restricted to the loading range where the papers report small modulation | gas and particle mean velocity within the stated measurement uncertainty; limits of the model written next to the result |
| V28.19 | SP5 | Evaporating spray with measured droplet statistics by size class at several stations (mean and fluctuating velocity, mass flux, diameters), through datums with histograms and size classes (D264) | measured: Sommerfeld and Qiu (1998), phase-Doppler measurements in a turbulent flow; **liquid, boundary conditions, data access and licence to be read first**; a liquid other than water needs a property set and a `DATA_SOURCES.md` entry | profiles within the measurement uncertainty or the difference written down with its probable cause; no tuning of model constants to the data (constants stay at the literature values, D201) |
| V28.20 | SP5 | Single droplet evaporation against measurements: droplet diameter squared in time, in still and convective hot air, for water | measured: Ranz and Marshall (1952), Hubbard, Denny and Mills (1975); model review: Sazhin (2006); exact d² law (Godsave 1953, Spalding 1953) as the limit | inside the scatter of the data for the range of the papers; the d² law's limit reproduced exactly (V28.7) |
| V28.21 | SP5 | Plug-flow reference of the spray-cooling duct (D240): droplets of one size in a uniform stream, momentum, heat and evaporation, ODE system integrated at 1e-10; compared with the CFD duct's mass-flow-weighted averages at the datums along the duct, and the outlet temperature against the adiabatic-saturation limit | own derivation from the cited laws (Ranz-Marshall, Abramzon-Sirignano, Spalding); exact states from the enthalpy balance | profiles along the duct within the band set by the plug-flow assumptions, stated before the run (the CFD has profiles, the ODE has not); the outlet within the error band of V28.4; the exact limit reached when water is in excess (V28.5, V28.11) |

When each check is run: V28.12 to V28.15 with SP1b; V28.1 to V28.4, V28.6,
V28.16 to V28.18 with SP2; V28.5, V28.7, V28.8, V28.19 to V28.21 with SP4 and
SP5; V28.9 to V28.11 with SP2s. A check that fails is reported with its
numbers and its probable cause; a band is never widened after the run to make
it pass (the rule of the verification grades, D199).

## 4. Staging (decided 2026-11-24: datum with SP1b, the rest per stage)

| Stage | Content | Starts when |
|---|---|---|
| SP1b | Datum object for gas and spray (D262): planes, boundaries, face zones, surface files, regions; gas reports (mass flow, area and mass-flow averages, enthalpy flow); liquid crossings and cumulative flows (D263); moment sums, diameters and histograms (D264); spray-wide numbers that need no coupling (D265); `datums.csv`, histograms, `result.json`, run record, tools and Python reader (D267); `series.json` manifest (D268); checks V28.12 to V28.15; gallery example, case, recipe with images | **Implemented 2026-11-24 (#265)** except the enthalpy flow (SP5) and V28.12; V28.13 to V28.15 pass ([gallery 265](../gallery/sprays.md#spray-datums)) |
| SP2, SP5, SP2s | Averaged fields (D266) with SP2; evaporation rate, evaporation length, vapour penetration with SP5; the checks of the table above at their stage | their own go |
| MurmelCraft M1 | Parcel viewer and time player (D268), in the MurmelCraft repository; can start as soon as the parcel file layout of SP1 is merged | go |
| MurmelCraft M2 | Charts, tables and statistics panel | go, after M1 |

## 5. Open points

- **Turbulence modulation by particles** (extra production and dissipation in
  k and epsilon from the parcel exchange) is in no decision of the
  specification. It matters for V28.18 and for dense sprays. Own decision when
  SP2 starts; until then the limit is written next to every comparison.
- **Which two-way measurement** (V28.18) and the exact conditions of V28.19
  are chosen after reading the papers, with the maintainer.
- **Cut-plane accuracy** of D262: the first version uses the cell value and
  gradient; a flux-conservative cut (using the face fluxes around the plane)
  is a later refinement if V28.13 shows an error that matters.
- **Statistical error of the time-averaged fields** in cells with few parcels
  (D254 (5), D266): a minimum-sample mask is the candidate; decided with the
  first gallery case.

## 6. References (cited from the literature record, not yet read)

- Csanady, G. T. (1963), Turbulent diffusion of heavy particles in the
  atmosphere, J. Atmos. Sci. 20.
- Gunn, R., Kinzer, G. D. (1949), The terminal velocity of fall for water
  droplets in stagnant air, J. Meteorology 6, 243-248.
- Hubbard, G. L., Denny, V. E., Mills, A. F. (1975), Droplet evaporation:
  effects of transients and variable properties, Int. J. Heat Mass Transfer 18.
- Kulick, J. D., Fessler, J. R., Eaton, J. K. (1994), Particle response and
  turbulence modification in fully developed channel flow, J. Fluid Mech. 277.
- Lau, T. C. W., Nathan, G. J. (2014), Influence of Stokes number on the
  velocity and concentration distributions in particle-laden jets, J. Fluid
  Mech. 757.
- Ranz, W. E., Marshall, W. R. (1952), Evaporation from drops, Chem. Eng.
  Prog. 48.
- Sazhin, S. S. (2006), Advanced models of fuel droplet heating and
  evaporation, Prog. Energy Combust. Sci. 32.
- Sommerfeld, M., Qiu, H.-H. (1998), spray evaporation in turbulent flow with
  phase-Doppler measurements, Int. J. Heat Fluid Flow 19 (title, liquid and
  data access to be confirmed when read).
- Snyder, W. H., Lumley, J. L. (1971) and Wells, M. R., Stock, D. E. (1983):
  as in the [assessment](../assessments/spray-module-prerequisites.md).
- Crowe et al. (2011), Godsave (1953), Spalding (1953), Abramzon and
  Sirignano (1989): as in the specification and the assessment.
