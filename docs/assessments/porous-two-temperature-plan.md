# Two-temperature porous module: plan (D270 to D279)

Plan, 2026-11-24. Requested by the maintainer after the honeycomb thermal
storage walkthrough of the same day (card 3: "two-temperature porous module,
plan first"; card 4: "both: correlations by default and tables from resolved
unit-channel runs"). Documentation only; **coding of each stage starts on the
maintainer's explicit go.** State read: `main` at 33de454, `docs/requirements.md`
(D18, D62, D85, D88 to D92, D99, D110, D112, D116, D118, D119, D150, D158,
D160, D201, L34, L48), the [radiation specification](../requirements/step29-radiation.md)
(R1a to R1g, D255 to D259, D269), the [MPI sizing](mpi-multi-rank-sizing.md)
and the crates `foehncraft-flow` (`case.rs`, `materials.rs`, `energy.rs`,
`cht.rs`, `frozen.rs`).

Clean room (D201): methods from papers and textbooks, correlations and exact
solutions as references. Results of other solvers (including the maintainer's
earlier two-dimensional porous model of the storage) are not used as reference
numbers without a decision of their own.

## 1. Answer in short

The target is a thermal storage built from ceramic honeycomb bricks (the
"Solarstone" brick: 150 × 150 × 150 mm, 11 × 11 square channels of 6.9 mm),
charged by hot gas from one end and discharged from the other, with the
moving thermocline in the solid as the main result and the heat capacity per
brick volume as the design goal.

Two ways to compute it, both in the plan:

1. **Resolved unit channel (needs radiation R1b).** One quarter of one
   channel with its share of the wall, cut by two symmetry planes, solved as
   conjugate heat transfer. The cut is valid: the brick repeats the channel in
   both cross directions, and for radiation a symmetry plane is a specular
   mirror (R1a), which is exactly the image of the neighbouring channels. This
   is the R1b honeycomb case (maintainer, card 2).
2. **Porous block in a large domain (this plan).** The whole brick stack is
   one porous region with two temperatures per cell, gas and solid. Its
   parameters come from correlations by default, or from tables that the
   resolved unit-channel runs of 1. write over Reynolds number and
   temperature.

What exists today and what is missing:

| Piece | Today | Missing |
|---|---|---|
| Region kind `porous` | parsed (D112), **refused** (REG001, `case.rs`) | momentum sink, superficial velocity, porosity |
| Directional resistance | material frames of orthotropic k exist (`materials.rs`, L48 axes) | a resistance tensor in those frames |
| Solid temperature in the porous cells | — | a second temperature in the energy matrix, exchange h·a |
| Anisotropic solid conduction | tensor k on solid regions (D116, gallery "orthotropic patch") | effective k of the honeycomb, radiative part k(T) |
| Long charge cycles | frozen-flow march (D118, #208), transient CHT | charge/discharge switching by `[[schedule]]`: **known bug**, the schedule and the re-read drop the energy boundary conditions (`case.rs`, `schedule.rs`), so a cycle with an inlet temperature fails validation |
| Large temperature range, mass-flow inlet | variable density (7c) | 7c refuses mass-flow inlets (ENE006); a storage is charged by a mass flow |
| Ceramic material data | user tables (D100) | no ceramic in the database; needs a datasheet or a source with a D201 entry |
| Unit-channel tables | study layer (D62), Nu and f per module (D110) | a pipeline that turns runs into a porous table; radiation in the unit channel (R1b) |
| Outputs | datums of the spray page (D262 to D265), ledger accounts (D119, #238) | thermocline position and thickness, stored heat, utilisation |
| Ranks, GPU | flow on ranks (W5), energy on ranks with W6 | the porous terms are cell-local; follow the energy equation |

Sizing (agent time): five stages of about 30 to 90 minutes of code each,
verification runs of minutes, the upscaling validation one to two hours of
runs; CI separately (section 9).

## 2. Brick numbers used throughout

For the 150 mm brick with 11 × 11 channels of width w = 6.9 mm:

| Quantity | Value | How |
|---|---|---|
| pitch p | 13.6 mm | 150 mm / 11 |
| porosity ε (open frontal area) | 0.256 | 121 · w² / (150 mm)² |
| hydraulic diameter d_h | 6.9 mm | square channel |
| specific surface a | 148 m²/m³ | 4 ε / d_h |
| brick length in channel diameters | 21.7 | 150 / 6.9 |
| Reynolds number in a channel | about 100 to 1 000 | air at 300 °C, 0.5 to 5 m/s in the channel: laminar |
| Graetz number x* = L / (d_h Re Pr), Re = 430 | about 0.07 | thermally developing over much of the brick |

The meaning of the "690" in the brick name is not known here (the maintainer
can confirm). The last row matters: with x* near 0.07 the fully developed
Nusselt number 2.98 under-predicts the mean heat transfer of one brick, and
every brick joint can restart the boundary layers (D273, Q1).

## 3. Porous flow (D270)

**Formulation.** Superficial (Darcy) velocity u as the unknown in the porous
cells, so that mass flux is continuous across the porous–clear interface;
volume-averaged momentum (Whitaker 1999; Nield & Bejan 2017, ch. 1 and 2;
Hsu & Cheng 1990):

  ∂(ρu)/∂t / ε + ∇·(ρ u u / ε²) = −∇p + ∇·τ_eff − (μ K⁻¹ + ½ ρ |u| F) · u + ρ g

with the permeability tensor K⁻¹ and the Forchheimer tensor F given by
principal values in a **material frame** (the axes of L48, the same input as
the orthotropic conductivity of today: `axes = [[..], [..]]`). The convective
and the Brinkman (viscous) terms are kept, because they matter at the
interfaces and are cheap; key `brinkman = true` (default).

**Along the channels.** From the laminar square-duct friction factor
f_D Re = 56.9 (Shah & London 1978) and u_channel = u / ε:

  1/K_∥ = (f_D Re / 2) / (ε d_h²) = 28.45 / (ε d_h²),  F_∥ from the apparent
  friction of the developing flow and the entrance and exit losses (D274).

**Across the channels** the walls are closed. The across resistance is a
ratio to the along value (`blocking = 1e4` default, a key), implicit in the
momentum matrix, so the velocity across is effectively zero without a
special constraint and any orientation of the brick axes works. The exact
alternative, a velocity projected onto the channel axis, is rejected for the
first version: it breaks the pressure–velocity coupling at the interfaces and
needs its own Rhie–Chow form.

**Discretisation.** The sink goes into the diagonal (Patankar) for SIMPLEC
and into the momentum block of the coupled solver (D176); the Rhie–Chow
face velocity uses the porous diagonal on both sides, and the porosity jump
at an interface is treated as a face property (harmonic as D99), so a
uniform flow through a block gives no pressure wiggles (verified, V1). The
Forchheimer term is linearised with |u| from the last iterate.

## 4. Two-temperature energy (D271)

Gas temperature T_g and solid temperature T_s per porous cell, not forced
equal (local thermal non-equilibrium; Schumann 1929; Kaviany 1995, ch. 3;
Nield & Bejan 2017, sec. 2.2.2):

  gas:   ε ρ_g ∂h_g/∂t + ∇·(ρ_g u h_g) = ∇·(k_g,eff ∇T_g) + h a (T_s − T_g)
  solid: (1 − ε) ρ_s ∂h_s/∂t = ∇·(K_s,eff ∇T_s) + h a (T_g − T_s) + q_s

with the conserved enthalpies of D92 and the T-correction form of D150 for
both phases.

**Monolithic, as D158.** The gas temperature of a porous cell is the ordinary
temperature of the energy domain (continuous with the clear fluid). The solid
temperature of each porous cell adds **one extra row** to the energy matrix,
coupled to the gas row of the same cell by h·a·V and to its porous neighbours
by the solid conduction. So the whole storage, gas, solid bricks, housing and
insulation is still **one temperature matrix** solved by the existing
BiCGStab plus aggregation AMG on the domain view. The exchange h·a is large
(stiff: h ≈ 15 to 30 W/(m² K), a = 148 m²/m³), and implicit coupling keeps
the time step free of it. The AMG strength is measured on the scaled matrix
(D158), so aggregates do not mix the two phases where the coupling is weak.

Rejected: a segregated solid equation with a point-implicit exchange. It
converges slowly when h·a is large and loses the exact conservation of the
exchange term within one iteration.

**Frozen flow (D118).** Charging lasts hours; the flow changes only when the
mass flow or the gas temperature changes much. The frozen-flow march of
#208 carries both temperatures with large time steps, and the flow is
refreshed when the density changes by more than a set fraction (key).

## 5. Effective conductivity of the solid and of the gas (D272)

**Solid, exact along, estimated across.** The walls of the brick are
continuous in all three directions:

- along the channels: K_s,∥ = (1 − ε) k_s, exact;
- across: the walls parallel to the heat flow form strips of fraction
  1 − √ε, the walls across add in series: K_s,⊥ ≈ k_s (1 − √ε) plus a small
  gas term (a parallel–series estimate; Kaviany 1995).

So for this brick (ε = 0.256) the solid conducts about 0.74 k_s along and
0.49 k_s across: anisotropic, but the walls do not make the across direction
an insulator. (A first summary of 2026-11-24 said "low across"; the numbers above
correct that.) Bricks stacked with gaps or misaligned
joints lower the along value; a joint conductance per brick length is a key
(default none).

**Radiation inside the channels.** At 400 °C and above, the radiation along a
channel carries heat from hot to cold wall sections. It is added to the
solid along the axis as a radiative conductivity

  k_rad,∥(T_s) = c_rad · 4 σ T_s³ d_h,

with the factor c_rad from the unit-channel tables (D273; it depends on the
wall emissivity and on L/d_h) and **zero until R1b** gives it. Across the
channels radiation goes through the walls only by conduction, so no
radiative part there. This keeps radiation out of the porous solve: no Monte
Carlo in the large domain.

**Gas.** ε k_g in all directions. Thermal dispersion along the channels
(Taylor–Aris) is small at these Péclet numbers compared with the h·a
exchange; an optional additive along-axis value (key, default off).

## 6. Parameters: correlations by default, tables from unit-channel runs (D273)

Decided by the maintainer on 2026-11-24 (card 4): both sources.

**Correlations (default)**, laminar square channels (Shah & London 1978):

- friction: f_D Re = 56.9 fully developed; apparent friction with the
  developing length (incremental pressure drop K(∞) about 1.43);
- heat: Nu_T = 2.98 fully developed at constant wall temperature; the mean
  Nusselt number of the thermally developing flow over one brick length as a
  function of x*, boundary layers restarting at each brick joint (decided,
  Q1);
- h = Nu k_g / d_h, a = 4 ε / d_h from the brick geometry keys (`channel`,
  `pitch`, `brick_length`), so the user gives geometry and gets K, F, h·a,
  K_s,eff.

**Tables from resolved runs.** The study layer (D62) runs the quarter channel
(section 1, way 1) over a grid of Re (5 to 8 values) and temperature (3 to 5
levels, with radiation once R1b exists) and writes one table file:
`f_app(Re)`, `Nu_m(Re, T)`, `c_rad(T)`, each with its mesh hash (D85), case
hash and solver version as provenance. The porous region reads it with
`parameters = { table = "solarstone.toml" }`; values in between are
interpolated (linear in log Re), outside the range the run stops with a
diagnostic (D78) instead of extrapolating. The resolved runs use the
D110 evaluation (bulk temperature, Nu and f from face fluxes). Periodic meshes
with solids are refused today; the quarter channel therefore uses one brick
length with an inlet plenum, which also gives the entrance effect that the
porous model needs.

## 7. Interfaces (D274)

- **Porous to clear gas** (plenums above and below the stack): superficial
  velocity and pressure continuous; optional contraction and expansion loss
  coefficients K_c, K_e at the face from the open area ε (Kays & London 1984),
  booked as a face pressure jump. The gas temperature is one field across the
  interface. The solid end faces exchange heat with the gas only through h·a
  inside the cell (default) or by radiation to the plenum walls once R1b is
  coupled to porous ends (later).
- **Porous to a solid region** (housing, insulation): the solid temperature
  of the porous cell conducts to the housing through the interface conductance
  of D160, with the porous solid share (1 − ε) of the face; the gas at the
  housing has a slip wall (the Darcy term dominates) and an adiabatic face.
- **Porous to porous** (two brick types): both temperatures continuous;
  resistance tensors as face properties (D99).

## 8. Outputs, ledger and the second law (D275, D276)

**Outputs (D275)**, as structured output and in the run record (D62), every
N iterations or time steps:

- stored heat Q_s = ∫ (1 − ε) ρ_s (h_s(T_s) − h_s(T_ref)) dV of the porous
  solid, and its rate;
- thermocline position (the axial position where T_s is the mean of the hot
  inlet and the cold initial temperature) and thickness (10 % to 90 % length),
  per porous region along its material axis;
- utilisation (stored heat over the heat of the brick volume at the inlet
  temperature);
- pressure drop over the stack and gas outlet temperature over the cycle;
- charge and discharge efficiency over a full cycle.

**Ledger (D276)**: a new storage account "porous solid" per region (D119
module rows, #238), the exchange pair gas–solid booked on both sides (net
zero to round-off), and the entropy production of the exchange
h a (T_g − T_s)² / (T_g T_s) ≥ 0 as a second-law check (D90).

## 9. Ranks and GPU (D277)

The porous terms are cell-local: the sink, the exchange and the extra solid
rows live with their cell, so the distribution of W3 to W5 holds without new
communication beyond one more field in the halo exchange (T_s). The extra
rows are owned by the rank that owns the cell; global ids of the solid rows
follow the cell global ids, so the reproducible mode (D93) is kept. Sums of
the outputs and the ledger use the exact sums of W2. GPU: the sink and the
exchange are pointwise kernels; the extended matrix uses the existing device
path of the energy equation. Ranks and GPU come in PM5, after the serial
stages (Q2).

## 10. Verification and validation (D278)

Exact or analytic, no other solver:

1. **V1 Darcy and Forchheimer**: one-dimensional flow through a block:
   pressure gradient exact to round-off; a block in a duct with clear
   plenums gives no wiggle at the interfaces.
2. **V2 directional block**: the brick axes turned by 30° in a duct; the
   flow follows the axis, the across velocity scales with 1/blocking.
3. **V3 exchange relaxation**: closed cell, no flow, T_g ≠ T_s at start: the
   difference decays as exp(−h a (1/C_g + 1/C_s) t), exact.
4. **V4 Schumann (1929)**: step change of the inlet gas temperature into a
   bed without conduction: the exact gas and solid temperatures in x and t;
   second-order convergence in space and time.
5. **V5 equilibrium limit**: h·a → large gives the one-temperature mixture
   solution, which the energy equation of today computes directly.
6. **V6 ledger**: stored heat equals the integral of the inflow minus outflow
   enthalpy to round-off; the entropy production non-negative.
7. **V7 upscaling consistency** (with the tables): one resolved brick (3 × 3
   channels, full length) against the porous brick with its tables: outlet
   temperature history and stored heat within a set tolerance (proposed 2 %).
8. **V8 measured data** (if available, see Q3): a storage charge with
   measured outlet temperature, with a D201 entry first.

Benchmark rule (geometry topics, 2026-11-21): the extra-row monolithic
solve is measured against the segregated point-implicit solve (iterations,
time per step) on V4 and the gallery case.

## 11. Gallery, cases and recipes (D18)

- **Gallery** (energy page, new topic "porous media"): the charging of a
  storage like the maintainer's image, **three-dimensional from the start**
  (Q4): a stack of several bricks with housing between two plenums,
  solid and gas temperature side by side at three times (Rainbow
  Desaturated, D261), the thermocline curve over time, the pressure field.
  No symbol-like patterns.
- **Cases**: `porous_darcy_block`, `porous_directional_block`,
  `porous_schumann`, `storage_charge_3d` (later `storage_cycle_3d` with
  charge and discharge).
- **Recipes**: "add a porous region with brick axes", "make porous tables
  from a unit channel" (the study-layer pipeline), "charge and discharge a
  thermal storage" (schedule, frozen flow, outputs).

## 12. Stages (D279)

| Stage | Content | Needs | Code (agent) |
|---|---|---|---|
| PM0 | Fix: `[[schedule]]` and re-read keep `boundary_energy` and `boundary_scalar`; mass-flow inlet under variable density (ENE006 part) | — | 30 to 45 min |
| PM1 | `kind = "porous"` lifted: superficial velocity, Darcy–Forchheimer tensors in material frames, blocking ratio, Rhie–Chow at porosity jumps, K_c/K_e; V1, V2; gallery: flow through a brick in a duct | — | 60 to 90 min |
| PM2 | Two temperatures: solid rows in the energy matrix, h·a, effective solid tensor, ledger account and entropy check; V3 to V6; benchmark against segregated | PM1 | 60 to 90 min |
| PM3 | Correlations (Shah–London, brick geometry keys), table reader, outputs of D275, frozen-flow charging; gallery storage charge, cases, recipes | PM2, PM0 | 45 to 60 min, runs of minutes |
| PM4 | Upscaling pipeline: study-layer quarter channel → table; V7; radiative factor c_rad when R1b is in | PM3; R1b for c_rad | 45 min code, 1 to 2 h runs |
| PM5 | Ranks and GPU for the porous terms (Q2: after PM1 to PM4) | W6 | 45 to 60 min |

**PM1 as built (#273).** The sink is implicit in the momentum diagonal and
the off-diagonal parts of the tensor are deferred; the Rhie–Chow factor
takes the smallest principal value, because the factor is one number per
cell. That is enough for axes along the mesh axes, and the build refuses
turned axes (POR002): a turned honeycomb stalls with a scalar factor. A
tensor (or projected) Rhie–Chow factor is the follow-up that lifts it; the
3D storage of PM3 has its bricks along the mesh axes and does not need it.
Channels across the main flow converge slowly. K_c/K_e of the brick ends
are in `honeycomb.length` (K(∞) = 1.43 over the brick length). V1 holds
with second order at the block faces (0.41 % on 100 cells); a filled duct
is exact.

## 13. Questions and answers

All five answered by the maintainer on 2026-11-24 (18:36 to 18:40 UTC).

| # | Question | Answer |
|---|---|---|
| Q1 | Default heat transfer coefficient | **Developing flow per brick** (as proposed): the mean Nusselt number of the thermally developing flow over one brick length, boundary layers restart at each brick joint (D273) |
| Q2 | Ranks and GPU | **Serial first** (one rank with threads) for PM1 to PM4; ranks and GPU in PM5 (D277) |
| Q3 | Reference data and ceramic properties | **Both**: open-literature values of a cordierite-type ceramic with a D201 entry so that coding is not blocked; the maintainer's brick datasheet and any measured charge data replace them when available (D278 V8) |
| Q4 | First gallery case | **Three-dimensional stack from the start**, not the proposed 2D case: several bricks with housing between two plenums; runs of about 1 to 2 hours (section 11) |
| Q5 | Start before radiation | **Start PM0 to PM3 now** (as proposed), radiative conductivity zero; PM4 with c_rad when R1b is in. This is the maintainer's go for coding PM0 to PM3 |

## 14. References

Hsu, C. T. & Cheng, P. 1990, Thermal dispersion in a porous medium, Int. J.
Heat Mass Transfer 33(8). Kaviany, M. 1995, Principles of Heat Transfer in
Porous Media, 2nd ed., Springer. Kays, W. M. & London, A. L. 1984, Compact
Heat Exchangers, 3rd ed., McGraw-Hill. Nield, D. A. & Bejan, A. 2017,
Convection in Porous Media, 5th ed., Springer. Schumann, T. E. W. 1929, Heat
transfer: a liquid flowing through a porous prism, J. Franklin Inst. 208.
Shah, R. K. & London, A. L. 1978, Laminar Flow Forced Convection in Ducts,
Academic Press. Whitaker, S. 1999, The Method of Volume Averaging, Kluwer.
