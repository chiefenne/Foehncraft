# Industrial steam furnace: development path

Plan, 2026-11-24. Requested by the maintainer: "For the solver I want later to
be able to calculate an industrial furnace for steam production. That requires
gas reacting flows, radiation, conjugate heat transfer, multi-phase boiling.
Probably more. Design a development path and ask me for decisions."
Documentation only; **coding of each stage starts on the maintainer's explicit
go.** State read: `main` at 9335cf1, `docs/requirements.md` (L1 to L3, L12, L14,
L15, L28, L31, L34 to L36; D17, D47, D88, D92 to D94, D102, D103, D107, D112,
D118, D119, D120, D194, D201, D238, D239), the [radiation specification](../requirements/step29-radiation.md)
(R1a to R1g, D255 to D259, D269), the [spray specification](../requirements/step28-sprays.md)
(SP1 to SP7, D240 to D254, D260), the [porous plan](porous-two-temperature-plan.md)
(PM0 to PM5, D270 to D279), the [MPI sizing](mpi-multi-rank-sizing.md) (W1 to
W9), the [graphics algorithms assessment](graphics-algorithms-for-modules.md)
and the multi-physics architecture plan of 2026-11-16 (D112 to D122).

Clean room (D201): methods from papers and textbooks, measurements and exact
solutions as references. Results of other solvers are not reference numbers.
Mechanism and property data enter through the importers of D103 with their
provenance and licence.

Decision numbers D280 to D285 are **provisional**; they are fixed by the
maintainer's answers (section 5) and checked against `main` before merge.
Stage ids use the prefix **SF** (steam furnace).

## 1. Answer in short

A fired steam generator is a chain of four coupled problems, and the solver
already has the middle of it:

1. **Gas side, burner and furnace chamber:** turbulent non-premixed reacting
   flow of natural gas and air, low-Mach variable density (7c, done), species
   and heat release (missing: L2), NOx as a post-step on the frozen flow
   (D118 exists).
2. **Radiation:** the dominant heat transfer to the water walls (70 to 90 %
   of the absorbed heat in the furnace chamber). Grey Monte Carlo kernel R1a
   is on `main` and its lead case is already a box furnace with cooled tubes
   (D255); missing: coupling R1b, the H₂O/CO₂ gas model R1f, scattering R1c
   for soot or particles.
3. **Walls and tube banks:** membrane water walls, superheater, reheater and
   economizer banks. Conjugate heat transfer with many fluids and solids is on
   `main` (D112, D194, thin walls and contact resistances in `cht.rs`); the
   two-temperature porous module (PM1 to PM5) gives the basis for tube banks
   as porous heat-exchanger zones. Missing: tube-bank correlations and a
   tube-side coupling.
4. **Water and steam side:** subcooled and saturated flow boiling in the
   tubes, natural or forced circulation, superheated steam. Missing entirely.
   IAPWS-95 water tables exist for the liquid; saturation, two-phase and steam
   states do not.

Recommended order (each a decision below): **gas side and radiation first**,
with the water walls held at the saturation temperature plus a fouling
resistance (a good approximation, because the inside of an evaporator tube is
within a few kelvin of T_sat); then a **one-dimensional tube-circuit model of
the water and steam side** coupled to the 3D walls, which is what the
industrial design of boilers needs; then **tube banks as porous zones**; the
**three-dimensional Eulerian two-fluid boiling** model (L3) comes after, for
local studies of a tube (departure from nucleate boiling, dry-out), not for
the whole furnace.

Agent time to the first reacting furnace with coupled radiation and water
walls (SF1 to SF5): **about 8 to 12 hours of code** in 8 to 10 pull requests,
plus validation runs of 1 to 2 hours each. To a full package boiler with
tube banks and the 1D steam side (SF1 to SF9): **about 16 to 25 hours of code**
plus runs. The 3D boiling model (SF10) is **not sized yet**; a first guess is
6 to 12 hours of code in several pull requests, with a high risk.

## 1a. High-end target (maintainer, 2026-11-24: "Go high end here!!!")

The path aims at an industrial-grade tool, not a demonstration:

- **Full package boiler scope** as the end point of this path (SF9): swirl
  burners, furnace chamber, membrane water walls, superheater, reheater and
  economizer banks, drum circuit with natural circulation, steam output,
  efficiency, tube metal temperatures, critical heat flux margin, NOx. Every
  earlier stage is a step towards that case, not a separate showcase.
- **Industrial-grade validation with measured data** at every level: laminar
  flames, turbulent jet and swirl flames (TNF workshop data), pilot-scale gas
  furnaces with in-flame temperature, species and wall heat-flux
  measurements, subcooled and saturated flow boiling, circulation and critical
  heat flux, tube-bank heat transfer and pressure drop, and full-scale boiler
  acceptance data where it is published. Clean room (D201) still applies:
  measurements and exact solutions only, no other solver's results.
  Each validation case gets an error band stated before the run and a verdict
  in `foehn compare` (L20).
- **GPU and MPI from the first stage**, as radiation did (D258, D259): every
  new kernel (species convection and diffusion, chemistry per cell, EDC,
  NOx, IF97 states, 1D circuits on rank 0 with exact sums, porous bank terms)
  is single-source for CPU and GPU with a CPU = GPU check, and every stage
  runs bitwise the same on 1 to 4 ranks. The MPI items W6 to W9 move onto
  the critical path (SF0), and GPU with ranks (v1.1) is needed by SF9.
- **Scale target:** a package boiler of 10 to 30 million cells on one
  multi-GPU node, and a utility boiler of 50 to 100 million cells on a
  cluster ("designed for" until measured, W9).

The added cost against the base plan: about 30 to 50 % more agent time per
stage for the GPU and rank paths and their checks (as R1a), plus the
validation runs, which are long (1 to 2 h each, several per stage). Revised
totals: SF1 to SF5 about 11 to 17 h of code; to SF9 about 22 to 35 h of code;
with the fuel scope of D281 (SF2h, SF11, SF13, SF15, SF16) about 35 to 55 h of
code, SF17 and SF10 not sized yet;
validation runs about 15 to 30 h spread over the stages. SF10 stays not sized
yet.

## 2. What the furnace asks for

| Zone | Physics | Scale | Solver part |
| --- | --- | --- | --- |
| Burners (swirl, staged air, flue-gas recirculation) | Turbulent non-premixed reacting flow, swirl, mixing | mm to 0.5 m | L2 reacting flows, RANS (SST, k-ε; curvature correction exists), LES optional (L5) |
| Furnace chamber | Hot gas, H₂O/CO₂ radiation, soot (weak for gas), NOx | 5 to 30 m | L15 radiation, L2 |
| Membrane water walls | Radiative and convective heat flux into fin-tube walls, fouling layer, boiling inside | tube 50 mm, fin 10 mm | CHT (D112), thin walls, 1D steam side |
| Convective pass (superheater, reheater, economizer) | Cross-flow over tube banks, flue-gas radiation in the cavities | thousands of tubes | porous zones with tube-bank correlations, or resolved sections |
| Steam side | Subcooled and saturated boiling, two-phase pressure drop, circulation, superheated steam, critical heat flux margin | tube circuits | 1D network (new), IAPWS-IF97 states |
| Drum, headers, attemperators | Mixing, spray attemperation | m | lumped volumes in the 1D network; sprays (L12) for attemperators later |
| Structure | Tube metal temperature, thermal stress, creep margin | mm | CHT, L36 thermoelastic (later) |
| Size of a run | 10 to 100 million cells for a utility boiler, 2 to 20 million for a package boiler | | MPI (W6 to W9), GPU with ranks (v1.1) |

## 3. What exists on `main`

| Need | Exists | Where |
| --- | --- | --- |
| Low-Mach variable density, p₀(t), buoyancy | yes | 7b1, 7c (#186, #195) |
| Energy, CHT with many fluids and solids, thin walls, contact resistance | yes | D150, D160, D194; `cht.rs` |
| Start values per region, absolute pressures | yes | D238, D239 |
| Turbulence RANS for swirl (SST, k-ε, curvature correction, QCR) | yes | stage 6, D236 |
| Bounded scalar convection (MLP limiter) for mass fractions | yes | D237 (#211) |
| Frozen-flow march (NOx and other post-steps) | yes | D118 (#208) |
| Ledger as data (element masses can be added as quantities) | yes | D119 (#238) |
| Material program, NASA-7/9 forms and importer, IAPWS-95 water liquid tables | yes, liquid only | D94, D102, D103; `materials.rs` |
| Monte Carlo radiation kernel, grey, multi-rank, GPU | yes (R1a, standalone) | #254, #255, #257, #259 |
| Lagrangian parcel layer (sprays, later soot or ash particles) | one-way (SP1) | #242 |
| Two-temperature porous zones | plan only | PM1 to PM5 |
| Multi-rank flow, energy, turbulence | yes, except solid regions, restarts, ζ-f | W5 (#234, #244); W6 next |
| Species transport, chemistry, heat release | **no** | L2 |
| Water and steam saturation, two-phase states (IAPWS-IF97) | **no** | D102 lists IAPWS; only liquid tables built |
| Boiling (any form) | **no** | L3; D123 F2 boiling curve for the spray wall film only |
| Tube-bank models | **no** | — |

## 4. What is missing, with dependencies

| Gap | Stage | Needs |
| --- | --- | --- |
| Species transport on the 7c density, mixture properties, element ledger | SF1 | 7c, D119 |
| Gas-phase chemistry program (subset of D120), stiff integrator with fixed elimination order (C14), global mechanisms, own equilibrium by Gibbs minimisation | SF2 | SF1 |
| Turbulence-chemistry interaction (eddy dissipation concept; tables later) | SF3 | SF2 |
| Radiation coupled to energy, WSGG for H₂O/CO₂, scattering | SF4 = R1b, R1f, R1c | R1a, D119, W5 |
| Furnace demonstration with water walls at T_sat + fouling resistance | SF5 | SF3, SF4 |
| Thermal and prompt NOx on the frozen flow | SF6 | SF5, D118 |
| IAPWS-IF97 saturation, two-phase and steam states | SF7a | D94 |
| 1D tube circuits (boiling, circulation, superheat), implicit coupling to 3D wall faces | SF7b, SF7c | SF7a, CHT interfaces |
| Tube banks as porous heat-exchanger zones with tube-side 1D | SF8 | PM1 to PM3, SF7b |
| Package boiler application on ranks | SF9 | SF5 to SF8, W6 to W9 |
| 3D Eulerian two-fluid flow boiling (wall boiling) | SF10 | L3 core, SF7a |
| Soot model and soot radiation | SF11 | SF3, R1c, R1f |
| Thermal stress of tubes | SF12 | L36 |
| Flamelet or progress-variable tables (non-adiabatic) | SF13 | SF2 |
| Resolved bubbles (VOF with phase change) | SF14 | L3 VOF |

## 5. Decisions for the maintainer (cards, one at a time)

| # | Question | Options | Recommendation | Decision |
| --- | --- | --- | --- | --- |
| Q1 (D280) | Priority: gas side or boiling first | (a) gas side and radiation first, walls at T_sat; (b) boiling first; (c) both in parallel threads | (a): radiation is the lead physics of the furnace, R1a's lead case is already a furnace (D255), and an evaporator wall at T_sat is a sound boundary condition until the steam side exists | **Answered 2026-11-24 (D280): (c) both in parallel threads.** The gas side (SF1 to SF6 with R1b, R1f) and the steam side (SF7a to SF7c) are built in two threads at the same time; they touch different crates and meet in SF5/SF7c (the walls switch from T_sat to the 1D circuit) |
| Q2 (D281) | Fuel and chemistry scope of the first version | (a) natural gas, eddy dissipation concept with global mechanisms (2-step, 4-step); (b) natural gas, flamelet tables first; (c) natural gas with hydrogen blends from the start; (d) also oil and solid fuels | (a): the chemistry program of D120 is needed anyway for L28 and L31; the eddy dissipation concept works with a stiff integrator that we can make reproducible (D93) and run on the GPU; tables (SF13) later. Hydrogen blends are a small step after (a) | **Answered 2026-11-24 (D281): the most advanced scope**, read as the union of the options: gaseous fuels (natural gas with hydrogen blends up to pure hydrogen, with differential diffusion), fuel oil (L12 sprays with multicomponent evaporation) and solid fuels (coal and biomass particles: drying, devolatilisation, char oxidation, ash). Turbulence–chemistry: eddy dissipation concept with finite-rate chemistry (SF3) **and** non-adiabatic flamelet / progress-variable tables (SF13) in the first version. Soot (SF11) becomes part of the first version, because oil and coal flames are sooting. Stages SF15 to SF17 added |
| Q3 (D282) | Water and steam side model | (a) 1D tube circuits with boiling correlations, coupled implicitly to the 3D walls; (b) 3D Eulerian two-fluid boiling for the tubes; (c) VOF with phase change | (a) first: this is what boiler design needs (circulation, pressure drop, critical heat flux margin, metal temperature); (b) later for local tube studies; (c) research only | **Answered 2026-11-24 (D282): all three.** (a) 1D tube circuits (SF7a to SF7c) carry the boiler; (b) 3D Eulerian two-fluid flow boiling (SF10) for local tube studies (departure from nucleate boiling, dry-out, void profiles); (c) VOF with phase change (SF14) for resolved bubbles and the closure data of (b) (bubble departure, wall heat-flux partition), kept as a research path off the critical path. Order inside the steam-side thread: SF7a → SF7b → SF7c → L3 core → SF10 and SF14 |
| Q4 (D283) | Tube banks of the convective pass | (a) porous heat-exchanger zones on the porous module, calibrated from resolved bank sections; (b) resolve every tube; (c) both as options | (a) by default with (b) for a calibration section: the PM4 upscaling pipeline is the pattern; thousands of resolved tubes cost 10 to 100 times the cells | **Answered 2026-11-24 (D283): (c) both as equal options.** Porous tube-bank zones (SF8) and fully resolved banks are both first-class choices per case, each validated against the bank correlations and measurements of section 7; the resolved banks also calibrate the porous zones (PM4 pattern). Resolved banks add SF8r below |
| Q5 (D284) | Extras in the first furnace version | (a) NOx post-step and fouling resistance per wall; (b) also soot; (c) also tube thermal stress; (d) all | (a): both are cheap and asked of every boiler study; gas flames are weakly sooting, so soot (SF11) can wait; stress needs L36 | open |
| Q6 (D285) | Validation data sources | (a) open literature only; (b) also data sets under a user agreement (for example OECD/NEA benchmarks); (c) also the maintainer's own plant data | (a): clean room (D201) and redistribution in CI are simplest; (b) case by case with a decision of its own | open |

## 6. Stages and sizing

Agent time is active work on code and tests, not elapsed time. Runs and CI
waits are separate. "Not sized yet" means there is no basis.

| Stage | Content | Needs | Code (agent) | Runs |
| --- | --- | --- | --- | --- |
| SF0 | Prerequisites (critical path, high-end target): MPI W6 (restart on ranks), W7 (run loop, monitors on ranks), W8, solid regions on ranks; R1a part 4 (#259); GPU with ranks (v1.1) before SF9 | W5 | as in the MPI plan (W6 1 to 2 h, W7 45 to 90 min); solid regions on ranks 45 to 90 min | CI |
| SF1 | **Species transport.** N mass fractions on the 7c density (N − 1 transported, the last by closure, or all N with a mass-conserving correction velocity), mixture c_p and h from NASA-7 (exists), μ and k by Wilke and Mathur rules, Fick diffusion with Lewis or Schmidt numbers, turbulent Schmidt number, MLP-bounded convection, element ledger rows (D119), inlet compositions in mole or mass fractions. Checks: exact 1D diffusion, mixing layer, sum of Y = 1 to round-off, element closure. Gallery: coaxial jet mixing (CH₄ into air), case, recipe | 7c, D119 | 60 to 90 min (1 PR) | minutes |
| SF2 | **Gas-phase chemistry.** The gas-phase slice of the chemistry program D120 (Arrhenius, third body, fall-off later), a stiff implicit integrator per cell with fixed elimination order (C14) for CPU = GPU and reproducibility (D93), global mechanisms (Westbrook–Dryer 2-step, Jones–Lindstedt 4-step), mechanism import from the YAML format of D103 (GRI-Mech 3.0 for checks), own chemical equilibrium by Gibbs minimisation. Checks: 0D constant-pressure reactor against the equilibrium and the adiabatic flame temperature, 1D laminar premixed flame speed of CH₄/air against measurements | SF1 | 90 to 150 min (2 PRs) | minutes |
| SF3 | **Turbulence–chemistry interaction.** Eddy dissipation concept (Magnussen 1981; Magnussen 2005 constants) with the SF2 integrator in the fine structures, eddy dissipation model as fallback; heat release into h_s; limits and diagnostics for stiff cells. Validation: piloted non-premixed jet flame (Sandia Flame D, TNF workshop data) and a bluff-body flame. Gallery: jet flame temperature, OH and CO fields | SF2 | 60 to 90 min (1 PR) | 1 to 2 h |
| SF4 | **Radiation for the furnace** = R1b (coupling, already specified), R1f (WSGG for H₂O/CO₂ with the coefficients of Bordbar et al. 2014, one grey gas per bundle), R1c (scattering, needed for soot and particles later) | R1a, D119, W5 | not sized yet in step 29; estimate 45 to 60 min per stage (3 PRs) | 1 h |
| SF5 | **Furnace demonstration.** A box furnace with one or more swirl burners, membrane water walls as CHT walls at T_sat (from the drum pressure, IF97 saturation line of SF7a or a constant) with a fouling resistance (existing `contact_resistance`), heat-flux maps on the walls, furnace exit gas temperature, ledger of heat release, radiation, wall heat and stack loss. Gallery, case, recipe; validation against a pilot-scale furnace with open measurements (Q6) | SF3, SF4 | 45 to 60 min | 1 to 2 h |
| SF6 | **NOx post-step** on the frozen flow (D118): thermal (extended Zeldovich with equilibrium O and OH), prompt (De Soete), assumed β-PDF in temperature; ppm at 3 % O₂ as output | SF5 | 45 to 60 min | minutes |
| SF7a | **Water and steam states.** IAPWS-IF97 regions 1, 2, 4 (and 3 if needed above 16.5 MPa) for h, s, ρ, c_p, T_sat(p), with the IAPWS formulations for μ and k; quality and two-phase mixture properties; checks against the IF97 verification tables | D94 | 45 to 60 min | — |
| SF7b | **1D tube circuits.** Graphs of tubes, headers and a drum: mass, momentum (friction, gravity, acceleration; Friedel two-phase multiplier; homogeneous or drift-flux slip), energy in h; heat transfer by Gnielinski (single phase), Chen or Gungor–Winterton (flow boiling), critical heat flux by the Groeneveld look-up table (2006) with margin as output; natural circulation by a fixed-point on the drum loop. Checks: exact single-phase cases, Bartolomei subcooled boiling (void and wall temperature), circulation ratio | SF7a | 2 to 3 h (2 PRs) | minutes |
| SF7c | **1D–3D coupling.** Each tube segment maps to its wall faces on the 3D mesh; the coupling is implicit through a conductance (h, T_fluid per face as a Robin condition, the interface pattern of D160 and D274), iterated to a common tolerance with the outer iterations; ledger pair between the 3D wall and the 1D circuit | SF7b, CHT | 60 to 90 min | 1 h |
| SF8 | **Tube banks as porous zones.** Porous zones (PM1 to PM3) with inline and staggered bank correlations (Žukauskas for heat transfer, Žukauskas or Gaddis–Gnielinski for pressure drop), anisotropic resistance in the bank axes (L48), gas-side h·a and the tube side from SF7b per cell (dual-cell coupling), cavity radiation through the gas (SF4); calibration from a resolved bank section by the PM4 pipeline. Gallery: superheater bank | PM1 to PM3, SF7b, SF4 | 2 to 3 h (2 PRs) | 1 to 2 h |
| SF8r | **Resolved tube banks as a first-class option.** Mesh generation for inline and staggered banks by script (Gmsh, thousands of tubes, periodic sections), tube-side 1D circuits on resolved tube walls through SF7c, local metal temperatures; checks against the same correlations and measurements as SF8; the porous zones read their tables from these runs | SF7c, SF8 | 90 to 120 min | 2 to 4 h (large meshes, ranks) |
| SF9 | **Package boiler application** on ranks: burner, furnace, water walls, superheater and economizer banks, drum circuit; steam output, efficiency, metal temperatures, critical heat flux margin, NOx. Application page, gallery images | SF5 to SF8, W6 to W9 | 60 to 90 min plus meshing | 2 to 4 h |
| SF10 | **3D flow boiling** (Eulerian two-fluid): two-fluid equations with phase pressure coupling, interphase drag, lift, wall lubrication, turbulent dispersion, bubble size by a one-group interfacial area equation, wall heat-flux partitioning (Kurul & Podowski 1990). Validation: Bartolomei, and a measured void profile from the open literature | L3 core, SF7a | not sized yet (first guess 6 to 12 h, several PRs) | hours |
| SF11 | Soot (two-equation, Leung–Lindstedt–Jones 1991) and soot absorption in the WSGG model | SF3, R1c, R1f | 60 to 90 min | 1 h |
| SF12 | Thermal stress of tubes and headers (one-way, L36) | L36 | not sized yet | — |
| SF13 | Non-adiabatic flamelet or progress-variable tables as a second model | SF2 | 90 to 120 min | 1 h |
| SF14 | VOF with phase change (resolved bubbles, research) | L3 VOF | not sized yet | — |
| SF2h | **Hydrogen blends.** Mixture-averaged diffusion (differential diffusion, Soret optional), detailed H₂/O₂ sub-mechanism, flashback-relevant laminar flame speeds; checks: H₂/air and CH₄/H₂ laminar flame speeds against measurements | SF2 | 60 to 90 min | minutes |
| SF15 | **Fuel oil.** On the spray module (SP2, SP4, SP5, SP6 breakup): multicomponent or continuous-thermodynamics evaporation of heavy fuel oil, pressure-swirl and air-blast atomiser injection, fuel vapour into the SF2/SF3 chemistry; validation: an open spray-flame data set with droplet and gas measurements | SP5, SP6, SF3 | 2 to 3 h (2 PRs) | 1 to 2 h |
| SF16 | **Solid fuels.** Lagrangian particles of coal and biomass: drying, devolatilisation (two competing rates, Kobayashi et al. 1977), char oxidation and gasification (kinetic–diffusion, Baum & Street 1971; intrinsic rates later), particle heat-up with radiation, swelling and shape for biomass; volatiles into the gas chemistry; particle absorption and scattering in the Monte Carlo radiation (R1c); ash particles; on the L28 particle–gas chemistry (D120). Validation: pilot-scale pulverised-coal furnace with open in-flame data | SP2, SF3, R1c, L28 | 4 to 6 h (3 PRs) | 2 to 4 h |
| SF17 | **Ash deposition, slagging and fouling as a computed layer** (particle sticking by viscosity or melt fraction, deposit growth, its thermal resistance and emissivity fed back to the CHT walls) | SF16, SF5 | not sized yet | hours |

Order and parallel work:

```
SF0 ─┐
SF1 → SF2 → SF3 ─┐
R1b → R1f ───────┴→ SF5 → SF6
SF2 → SF2h ; SF3 → SF13, SF11
SP5, SP6 → SF15 ; L28, R1c → SF16 → SF17
SF7a → SF7b → SF7c ─┐
PM1..PM3 ───────────┴→ SF8 → SF9
L3 core → SF10        (later)
```

SF1 to SF3 (gas) and SF7a to SF7b (steam side) touch different crates and can
run in two threads at the same time. R1b and R1f are already in the radiation
specification and follow R1a part 4.

GPU: SF1 to SF3 kernels are cell-local (chemistry per cell, no reductions),
so they go to the GPU like the existing cell kernels; the integrator must keep
a fixed elimination order to give CPU = GPU. GPU with ranks is v1.1.

## 7. Validation list (open literature or exact)

| Physics | Case | Kind |
| --- | --- | --- |
| Species | 1D diffusion of a binary mixture, sum of Y, element closure | exact |
| Chemistry | adiabatic flame temperature and equilibrium composition of CH₄/air | exact (own Gibbs minimisation) |
| Chemistry | laminar flame speed of CH₄/air at 1 atm, φ = 0.7 to 1.4 | measurements in the open literature |
| Turbulent flame | Sandia piloted flame D (Barlow & Frank 1998), TNF workshop data | measurements, free with citation |
| Furnace | pilot-scale natural gas furnace with in-flame and wall heat-flux measurements from the open literature | measurements (Q6) |
| Radiation | slab and cube (R1a cases), WSGG against line-by-line results published for H₂O/CO₂ slabs | exact or published |
| Steam states | IAPWS-IF97 verification tables | exact (standard) |
| Flow boiling | Bartolomei et al. (1982) subcooled boiling in a heated tube | measurements |
| Critical heat flux | Groeneveld et al. (2007) look-up table | correlation from measurements |
| Tube banks | Žukauskas (1972) heat transfer and pressure drop of tube banks | correlations from measurements |

## 8. Flaws of industry practice to avoid (for the challenger)

- **Order-dependent chemistry tables.** In-situ adaptive tabulation answers
  depend on the order of the queries, which breaks reproducibility (D93) and
  differs between ranks. If tabulation is used, it is a precomputed table
  (SF13), not built during the run.
- **Mass fractions clipped to [0, 1] after the solve**, so that elements do not
  close. Here: bounded convection (MLP) and the element ledger as a check.
- **Radiation called rarely with a frozen source**, so that the wall heat flux
  carries an unknown lag. Here: the linearised source of R1b and the
  statistical error as a health number.
- **Explicit, lagged 1D–3D coupling** that oscillates or needs heavy
  relaxation. Here: implicit conductance coupling (SF7c).
- **Isotropic porous tube banks**: a bank resists cross-flow very differently
  from flow along the tubes. Here: tensors in the bank axes (L48, D270).
- **Wall boiling that depends on the near-wall cell size.** Here: the wall
  partitioning reads the wall-layer kernel at a fixed y⁺ (as the wall
  functions of D173) and is checked on two meshes.

## 9. Risks

- Stiff chemistry in the eddy dissipation concept can stall the steady solver;
  mitigation: chemistry time-step limits, pseudo-time per cell, the coupled
  solver's step control (L45).
- Monte Carlo noise in the wall heat flux of a large furnace; mitigation:
  running averages over calls (R1b), null-collision tracking (R1e).
- The 3D two-fluid boiling model is a large new solver path (phase-coupled
  pressure) with many closure models; it is kept out of the critical path.
- Furnace geometry: an open pilot-scale geometry with measurements is needed
  (Q6); otherwise the gallery boiler is our own design with no measured
  reference.

## 10. References

- B. F. Magnussen (1981), eddy dissipation concept, 19th AIAA Aerospace
  Sciences Meeting; B. F. Magnussen (2005), the EDC with revised constants,
  ECCOMAS thematic conference.
- C. K. Westbrook, F. L. Dryer (1981), simplified reaction mechanisms for
  hydrocarbon fuels.
- W. P. Jones, R. P. Lindstedt (1988), global reaction schemes for
  hydrocarbons.
- G. P. Smith et al., GRI-Mech 3.0 (1999), free mechanism with citation.
- R. S. Barlow, J. H. Frank (1998), piloted CH₄/air flames C, D, E and F,
  TNF workshop.
- M. H. Bordbar, G. Węcel, T. Hyppänen (2014), a line-by-line based weighted
  sum of grey gases model for H₂O/CO₂ mixtures.
- G. A. De Soete (1975), prompt NO formation; Ya. B. Zeldovich (1946), thermal
  NO.
- W. Wagner et al. (2000), IAPWS-IF97.
- J. C. Chen (1966), flow boiling correlation; K. E. Gungor, R. H. S. Winterton
  (1986).
- L. Friedel (1979), two-phase pressure drop.
- D. C. Groeneveld et al. (2007), the 2006 CHF look-up table, Nucl. Eng. Des.
- V. V. Bartolomei, S. G. Bashkirov, V. M. Borisov (1982), subcooled boiling
  in a heated tube, Thermal Engineering.
- N. Kurul, M. Z. Podowski (1990), wall heat-flux partitioning.
- A. Žukauskas (1972), heat transfer from tubes in cross flow, Adv. Heat
  Transfer.
- K. M. Leung, R. P. Lindstedt, W. P. Jones (1991), simplified soot model.
- V. Gnielinski (1976), single-phase heat transfer in tubes.
