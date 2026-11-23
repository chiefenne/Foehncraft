# Stage 7 plan: energy, buoyancy, conjugate heat transfer (v0.7)

Planning result, 2026-11-20; kept up to date in §7. Moved into the
repository on 2026-11-21 (D200). Based on the stage 7 spec
`docs/requirements/stage7-energy-buoyancy-cht.md` (D150–D164, PR #60,
1ac75e6), step 21 (D172–D175), the stage 6 status and main at 9b0544b.

## 1. Status audit: what stage 7 needs against main

| Area | Spec | On main today | Gap |
|---|---|---|---|
| Energy equation | D150: h_s on fluid + solid, δT form, residual in h, D92 Newton T(h) | Nothing. Temperature exists only as a **passive scalar** (`scalar.rs`, D51/D116) with constant Γ | Whole energy module (kernels + flow), restart of h and T |
| Thermal wall law, Pr_t | D135 (stage 6) | **Done** for passive scalars (#85–#87): Kader, Jayatilleke, Pr_t 0.85 and Kays–Crawford; heated pipe and channel DNS validated | Reuse; repeat channel DNS with the energy equation as a regression |
| Streamwise-periodic heat | D109–D111 | Heat-flux mode of a scalar periodic pair (`PeriodicMode::HeatFlux`) | **Wall-temperature mode** missing (V7.14 Nu 3.657) |
| Density as a field | F1 | **Still one number** (`FluxInputs::rho`, `MomentumInputs::rho` in `foehncraft-kernels/src/flow.rs`); `materials.rs` refuses a varying density | Cell and face ρ, μ, k, c_p fields in flux, momentum and scalar kernels, CPU and GPU |
| Material database | D94–D103 | Architecture, forms (constant, polynomial, Sutherland, power law, tables with Steffen, ideal gas), `[materials]`, `[regions]`. Data: **ρ and μ only**, 5 fluids (air, water, engine oil, ethylene glycol, glycerol) | c_p and k for all, NASA-7/9, IAPWS for water, solids with orthotropic k, glycol–water mixtures, importers of D103, MAT006–MAT009 |
| Gravity, buoyancy | D153, D154 | **Nothing**: only a constant `[fluid] body_force` vector; the only `gravity` in the code is the pump head of `rotor.rs` | `[gravity]`, `[buoyancy]`, reduced pressure p′, balanced-force face flux, wall pressure (F2, F3) |
| Coupled solver | not in the spec (step 22 came later) | The coupled p–u algorithm (`coupled.rs`, D176–D178) now runs beside SIMPLEC | **New finding**: the balanced force and the variable density must go into both algorithms |
| Buoyancy in turbulence | D155: stage 6 leaves a `buoyancy_production` hook | **Hook is not there** (no S_u/S_p hook in k, ε, ω, ζ) | Add the hook in stage 7 (7e1) |
| Thermal BCs | D156 | Scalar BCs: value or flux only | temperature, heat flux, heat flow, convection, radiation, thin wall, backflow T, THM/BUO codes |
| Regions, CHT | D112, D157, D158 | `[regions]` parse; solids refused (REG001); flow solves on all cells (F4) | Domain views (flow on fluid cells, energy on all), interface faces, contact resistance, AMG over views |
| Non-conformal interfaces | D58 (moved into stage 7) | Imprinting exists **only for surfaces of revolution** (`imprint.rs`, sliding rotors, D186) | General static imprinting of planar/arbitrary interfaces, CHT002 |
| Partitioning | D159 | METIS one constraint, no weights (F5) | Two constraints, per-view halo plans |
| Ledger | D119, D161 | Fixed structs `MassBalance`, `EnergyBalance` (kinetic only) in `balance.rs` (F6, was to come with stage 6) | Ledger as data, thermal part, kinetic → thermal exchange, D90 parts 1–3 |
| Verification data | D162 | `DATA_SOURCES.md` has no stage 7 entry | Entries and licence checks for every V7 reference |
| Version | — | Workspace 0.5.0; **v0.6 not yet tagged** | — |

Summary: stage 7 starts from zero physics but on a good base (generic
scalar transport, property architecture, thermal wall law, periodic heat
flux, imprinting machinery). Two stage 6 prerequisites (ledger as data, the
buoyancy hook) were not delivered and move into stage 7.

## 2. Proposed PR plan (small steps)

The six PRs of D164 are cut into 17 smaller ones. Each PR: tests,
section in `docs/verification/stage7.md`, changelog fragment
`changelog.d/<PR>-<name>.md`, README "What works today", Python keys (D61),
GPU kernels with CPU = GPU check (see Q4). Physics PRs get a gallery entry
with 1–3 field images (no four-fold or symbol-like patterns), case files and
a recipe on the new topic pages `docs/{gallery,recipes,cases}/energy.md`
(and `cht.md` from 7d). Core PRs have no gallery entry.

### Lane 0: groundwork (core, no physics change)

| PR | Content | Check | Gallery / recipe |
|---|---|---|---|
| 7.0a | Ledger as data (D119): terms as named rows per region; mass and kinetic energy moved onto it | V7.25: every catalogue case bitwise unchanged, ledger numbers identical | none (core) |
| 7.0b | ρ, μ as cell and face fields in flux, momentum, scalar kernels (F1), CPU and GPU; constant fields from today's properties | V7.25 bitwise; GPU parity on 3 cases | none (core) |
| 7.0c | Thermal material data: c_p, k for the 5 fluids, NASA-7/9 forms + importer, water (Q5), solids starter set (Q6), orthotropic k form, glycol–water mixture, MAT006–MAT009; `foehn materials` plot of each property | unit tests against the sources' tabulated points; `DATA_SOURCES.md` entries | properties page: property curves of air, water, oil, aluminium (plots, no field) + recipe "define a solid material" |

### Lane A: energy on the fluid

| PR | Content | Validation | Gallery / recipe |
|---|---|---|---|
| 7a1 | Energy module, fluid only, constant ρ (a per-cell field, refreshed from T): h_s with δT form, T(h) Newton, conduction with k(T), c_p(T), μ(T) by per-cell refresh of the face viscosity each outer iteration, SGDH with the stage-6 Pr_t, thermal BCs temperature / heat flux / heat flow / convection, viscous heating, `[energy]` keys, ENE and THM codes, thermal ledger | V7.1 MMS on 4 mesh types, V7.2 Couette Br = 1, V7.3 adiabatic stirred box (ledger + entropy), channel DNS Nu as regression vs #86 | heated Couette T profile; lid cavity heating up (T field) |
| 7a2 | Wall-temperature periodic mode (D109 (2)), energy on periodic pairs | V7.14 laminar pipe Nu 4.364 and 3.657, order 1.8–2.2 | heated periodic pipe, T cross-section |
| 7a3 | Radiation to ambient, combined convection–radiation, thin wall, backflow T, quartic wall balance | V7.11 slab, T_w to 1e-12 in ≤ 10 iterations | recipe "outer surface losing heat by convection and radiation" (+ small T image) |

| 7t1 | Thin walls 1D (D189, added 2026-11-20): layers, heat capacity, sources, contact resistance on wall boundaries, TWL codes | V7.26 multi-layer steady, V7.27 transient slab | duct wall with insulation, T across the layers; recipe "insulate a wall without meshing it" |

### Lane B: buoyancy (needs 7a1)

| PR | Content | Validation | Gallery / recipe |
|---|---|---|---|
| 7b1 | `[gravity]`, Boussinesq, reduced pressure p′ and `pressure_kind`, balanced-force face flux and cell force (D154) in **SIMPLEC and coupled**, wall pressure from the reduced jump, BUO codes | V7.7 (a) stratified box at rest, max \|u\| ≤ 1e-10 √(gL), oblique g, 4 mesh types perturbed | recipe "switch on gravity"; static-pressure field of the box |
| 7b2 | Laminar natural convection benchmarks (moved to v1.1 on 2026-11-23) | V7.15 square cavity Ra 1e3–1e6 (de Vahl Davis, Le Quéré), V7.12 Ostrach plate, V7.9 Aung–Worku mixed convection | cavity T and velocity at Ra 1e6; plate boundary layer |
| 7b3 | Transient buoyancy | V7.13 transient MMS (BDF2, CN), V7.8 Rayleigh–Bénard onset Ra_c 1707.76, V7.16 cubic cavity (Tric et al.) | 3D cubic cavity iso-surfaces of T; Rayleigh–Bénard rolls (2D, rolls in a row, not cells in a grid) |
| 7c | Low-Mach variable density ρ(T, p₀), closed volume p₀(t), pressure work options | V7.7 (b) gas column, (c) density step; V7.17 large-ΔT cavity (Vierendeels et al.), mass to round-off | large-ΔT cavity density field |

### Lane C: conjugate heat transfer (needs 7a1; can run beside lane B)

| PR | Content | Validation | Gallery / recipe |
|---|---|---|---|
| 7d1 | Domain views (flow on fluid cells, energy on all, D112), REG001 lifted for solids, region interface faces with series conductance, one-sided gradients, AMG over the energy view, CHT001/CHT003 | V7.4 composite wall (ratios 1–1e4) exact to 1e-9, V7.5 composite cylinder, V7.6 contact of two semi-infinite bodies | composite cylinder T field; recipe "add a solid region" |
| 7d2 | Contact resistance, thin walls between regions, orthotropic k (tensor), `[[energy.sources]]`, CHT004 | V7.4 with resistance, V7.10 thick-walled periodic channel Nu 8.235 | **CHT showcase**: heat sink with 5 fins in a channel (T in solid and fluid, velocity) |
| 7d3 | General static non-conformal interfaces by imprinting (D58 beyond surfaces of revolution), CHT002 | D58 cases: Poiseuille and T-junction across a non-conformal plane, V7.4 on non-matching layers | fluid and solid meshed separately (mesh cut + T field) |
| 7d4 | METIS with two constraints, per-view halo plans (D159, threads only; MPI stays stage 8) | V7.24: views, halos and matrix identical in global order for 1, 2, 4 parts | none (core) |

| 7t2 | Shell conduction and baffles (D189), surface mesh W1 and wall-face halo W2 moved from L12, film slots on walls, thin walls and interfaces (D190) | V7.28 fin, V7.29 thin wall vs resolved CHT, V7.30 baffle exchanger, V7.31 ledger box, V7.32 film slots | plate heat exchanger plate T field; recipe "model a sheet-metal part as a shell" |

### Lane D: turbulent heat transfer (needs 7a1, 7b1)

| PR | Content | Validation | Gallery / recipe |
|---|---|---|---|
| 7e1 | `buoyancy_production` hook in k, ε, ω, ζ/φ (stage 6 debt), G_b with SGDH, C_ε3 (Henkes), ω source | V7.25 bitwise without gravity; V7.18 Tsuji–Nagano plate with SGDH (reported, no band) | turbulent vertical plate T field |
| 7e2 | GGDH tensor heat flux (C_θ 0.3), G_b with GGDH | V7.18 ±10 %, V7.19 Ampofo–Karayiannis cavity ±15 %, V7.21 turbulent periodic pipe with energy vs Gnielinski | turbulent cavity T and velocity |
| 7e3 | Heated backward-facing step and wall-function variants with variable properties | V7.20 Vogel–Eaton St peak and reattachment | heated step St and T field |
| 7e4 | D172 Pr_t closures (`kays-crawford-hp`, `wall-profile`), D173 variable-property wall treatment, semi-local y⁺ | H-0 unit, H-1 and H-2 heated pipe at high Pr (needs TU Graz numbers, Q7) | high-Pr pipe T⁺ profiles |

### Close

| PR | Content | Check |
|---|---|---|
| 7f | Dual-number Jacobians of all new kernels, ⟨Av,w⟩ test, Colab notebook, restart bit-identity in V7.6, verification manual, release v0.7 | V7.22, V7.23 |

Order and parallel work: 7.0a → 7.0b → 7.0c → 7a1 are sequential (each
changes the kernels the next one uses). After 7a1 the lanes B and C can run
in parallel; D after 7b1. 7a2, 7a3 fit beside B/C. 19 PRs (with 7t1 and 7t2 of D189/D190).

## 3. Questions to the maintainer (2026-11-20; go given, later decisions in §6 and §7)

Q1. **Tag v0.6 first?** Main is still 0.5.0. Recommend: tag v0.6 now (or
after the ζ-f₀ sink decision), so v0.7 is a clean stage-7 diff.

Q2. **Stage 6 debt into stage 7?** The ledger as data and the buoyancy hook
were not delivered. Recommend: yes, as 7.0a and inside 7e1.

Q3. **Finer cut**: 17 small PRs instead of the six of D164 (record as an
amendment of D164 in a docs-only PR). Recommend: yes.

Q4. **GPU per PR or at the end (7f)?** Recommend: each PR brings its GPU
kernel and a CPU = GPU check (the PTX guard test caught missed entries in
6f); 7f keeps only the Colab notebook and the release.

Q5. **Water properties**: IAPWS-95 is the reference but costly per cell.
Recommend: IAPWS-95 (and IAPWS 2008/2011 for μ, k) in the database code,
evaluated into Steffen tables at the operating pressure when the case is
resolved (D96); the kernel only reads tables.

Q6. **Solid starter set.** Recommend 8: aluminium (pure and 6061),
copper, carbon steel, stainless 304, titanium, fused silica glass, silicon,
and one polymer (PTFE), with data from open NIST/published sources only. Fluids add
glycol–water at arbitrary fraction. Maintainer to name others.

Q7. **Experimental and DNS data** (Tsuji–Nagano, Ampofo–Karayiannis,
Vogel–Eaton, Redjem-Saad, TU Graz tables): as for the channel DNS, only
published numbers in the repo, tables fetched by a script into a git-ignored
folder. Recommend: same rule; maintainer to supply what is not openly
downloadable (TU Graz for 7e4).

Q8. **Turbulent model coverage**: which models must meet the V7.18–V7.21
bands? Recommend: ζ-f₀ and EB k-ε with GGDH must pass; SST and k-ε are
reported without a band.

Q9. **Start order**: buoyancy lane first, CHT lane first, or both in
parallel after 7a1? Recommend: both in parallel; if only one at a
time, CHT first, because it is the industrial need (heat sinks,
electronics, the applications track).

Q10. **Confirm spec defaults**: reduced pressure at pressure boundaries
when gravity is set (`pressure_kind = "static"` as option), C_θ = 0.3,
transformed ω source. Recommend: confirm; the validation cases may change
C_θ and the ω source.

Q11. **Go**: coding starts only on the maintainer's go. Recommend: go for
7.0a now (core, no physics change, no conflict with open PR #131 or the
ζ-f₀ work).


## 4. Added 2026-11-20: thin walls and film heat exchange

Recorded in PR #132 (D189, D190, O24 with TW1–TW6, checks V7.26–V7.32, PR 7t of D164). 7t1 after 7a1, before full CHT; 7t2 after 7d1.

## 5. Iteration budget (maintainer rule 2026-11-20)

Every validation run, steady ones first: look at the run after 2000-3000 iterations (what goes wrong, which settings to change); never let a run go to 10000+ iterations. Case files of V7.x set max_iterations accordingly.

## 6. Decided 2026-11-20 (maintainer, after the CHT question)

- Fluid properties temperature-dependent from 7a1 on: c_p(T), k(T), μ(T), each
  evaluated per cell after every temperature update (μ through the face
  viscosity buffer of 7.0b). Forms on the GPU: polynomial, table, power law,
  Sutherland, NASA-7/9; expressions in T stay on the CPU (no device-side
  program).
- ρ stays constant in 7a1, held as a per-cell field that can be refreshed
  from T; the Boussinesq ρ(T) of 7b reads the same field.
- Variable density with continuity and pressure correction (low-Mach, p₀(t),
  closed-domain thermodynamic pressure) stays in 7c.
- Shortest path to a first runnable CHT case: 7.0c → 7a1 → 7d1, laminar steady
  forced convection through a solid plate between two fluid streams. Thin
  walls, buoyancy, contact resistance and imprinting are not needed for it.

## 7. Status 2026-11-21

- 7a1 merged as #154 (energy on the fluid; walls through `[boundary.<name>.energy]`
  keys, THM001–005, THM009–011, ENE006; gallery and recipe 154). Moved on:
  V7.3, THM006–008, thermal ledger in `result.json`, device-side property refresh.
- 7d1 as #156 (first CHT): fluid view plus energy on all cells, series
  conductance, `interface-<region>` walls, generator `zones`. V7.4 exact; the
  first runnable case is `cht_heated_wall.toml` (heated solid wall under a
  laminar channel flow) instead of a plate between two streams. Still open
  from the 7d1 row: one-sided gradients, V7.5, V7.6, AMG over the energy view,
  CHT001/CHT003; steady laminar only.
- D204 (cycle-averaged CHT for periodic flows) recorded as #157; the
  maintainer moved turbulent and transient CHT before buoyancy (7b1).
- 7d5 as #158 (turbulent CHT: wall-function conductance at the interface;
  cases `cht_turbulent_channel`, `cht_heat_sink_turbulent`).
- 7d6 as #159 (transient CHT: one energy matrix with fluid and solid time
  terms; V7.6 exact, V7.37).
- 7d2 as #160 (contact resistance, thin walls, orthotropic solids, one-sided
  interface gradients, region sources).
- Open finding: a uniform start velocity into the front face of a solid
  makes the first SIMPLEC iterations fail; such cases start from rest
  ([pitfalls](../pitfalls.md#steady-runs)).
- Next in the plan: 7b1 (buoyancy), then lane D; v0.6 is tagged.
