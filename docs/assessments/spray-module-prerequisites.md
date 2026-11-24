# Lagrangian spray module (L12): what is missing before the first line of code

Assessment, 2026-11-24. Requested by the maintainer: "add the Lagrangian spray
module as soon as possible; which stages are missing so we can start?"
Discussion only, no code. State read: `main` at b431abe (v1.0.0 plus #224 to
#227), `docs/requirements.md` (L12, D47, D93, D107, D108, D117, D121 to D124,
D190), the wall film decisions in
[requirements](../requirements.md), the [MPI sizing](mpi-multi-rank-sizing.md) and the crates.

## 1. Answer in short

Nothing blocks a first, useful spray module. What exists is enough for
**non-evaporating liquid or solid parcels with drag, gravity, two-way momentum
coupling and turbulent dispersion on one rank** (spray cooling by sensible
heat only, scrubber-type and dust-type flows, painting without film). What is
missing for the next steps is, in this order:

1. A **written L12 v1 specification** (one page of decisions, D240 onward).
   Today L12 is one table row plus D47, D107, D117, D121 to D124; there is no
   model choice for drag, dispersion, injectors, output or restart of parcels.
2. The **particle layer** itself (new crate, structure-of-arrays, injection,
   tracking driver on `Locator`, deposit without atomics, output, restart).
3. **Multi-rank migration** (needs MPI W3 to W5; independent of 2).
4. A **minimal species transport (a slice of L2)** and a **liquid property
   data set**: both needed for evaporation, neither needed for stages S1 to S3.
5. The **surface mesh of the wall faces (D124 W1, W2)**: needed for the wall
   film only, not for stick or rebound.

Evaporation (with two-way energy coupling) is where the spray module becomes
what the maintainer means by "industry-standard"; it is gated by item 4.
Breakup is not gated by anything beyond S1 and S2.

## 2. What is on `main` and what it gives the spray module

| Capability | State on main | Use for sprays |
|---|---|---|
| Face walk through polyhedra | `crates/foehncraft-mesh/src/locate.rs`: `Locator::track`, `locate`, tetrahedral-decomposition walk (Macpherson, Nordin, Weller 2009), reports the exit face; point-location test on all mesh types | The tracking kernel. Reports the boundary face where a track leaves, so stick and rebound need no ray casting. Not yet used by the flow crate (tests only) |
| Counter-based RNG with stream word | `rng::keyed([stream, id, step, event])`; streams `SPRAY_INJECTION`, `SPRAY_BREAKUP`, `SPRAY_DISPERSION`, `SPRAY_WALL`, `SPRAY_COLLISION` already reserved (D121) | Done. Parcel ids from (injector, step, index) are not implemented (D117) |
| Transient solver, SIMPLEC and coupled p-u | merged | Time loop for a transient Euler-Lagrange run; the parcel step sits in the time step between the flow solve and the sources of the next solve |
| Energy equation, buoyancy and gravity | merged (7a1, 7b1) | Gravity vector for settling; gas temperature at the parcel |
| Low-Mach variable density, closed volume p0(t) | merged (7c, D151, D92) | Enough for sprays at atmospheric or elevated pressure and for evaporation into a closed box. Fully compressible flow (shocks, acoustics) is **not** needed for L12 and not present; sprays into high-pressure vessels work with the thermodynamic pressure |
| Turbulence k-eps, k-omega SST, zeta-f, EB variants | merged | k and epsilon (or omega) at the parcel for dispersion (Gosman and Ioannides 1983 form); no new turbulence work |
| Rotating frames: MRF, sliding interfaces | merged | Parcels in rotating zones need frame transformation and interface crossing: **not** in the first stages (see §6) |
| Scalar transport with face and cell fields, MLP limiter | merged | Reusable path for the vapour mass fraction and for smoothed parcel sources (D117 (2)) |
| Material database | forms exist (constant, polynomial, Sutherland, tables, ideal gas); rho, mu, c_p, k for five fluids | Gas properties at the parcel. **No** liquid data (vapour pressure, latent heat, surface tension, liquid c_p, liquid viscosity) and **no** mixture rule |
| BVH over wall faces | `query/bvh.rs`: closest point only (D108 part 1) | Enough for S1 to S3. Ray casting is D108 part 2 and is wanted with L15 and the film |
| Output | VTKHDF with one `PolyData` block per boundary (D75, D84) | Same format gets a `PolyData` block of parcels (points, diameter, velocity, temperature, parcel number) |
| Restart | single process, by global id; mapping to other rank counts is MPI W6 | Parcels are written by parcel id; the format needs a parcel block |
| Ledger | fixed structs `MassBalance` and `EnergyBalance` in `balance.rs`; the ledger as data (D119) is not done | The spray needs its own named exchange pair (gas momentum, mass, enthalpy against the parcels). With fixed structs this is one more struct; with D119 it would be rows. See decision Q4 |
| MPI | W1 (#200) and W2 (#219) merged; W3 to W5 being built; `decompose` refuses periodic meshes | Migration of parcels between ranks needs the processor faces and the halo of W3. Parcels on a rank must be reachable from the owner rank's cell list: the halo depth (D38) bounds the distance a parcel can move in one step |
| Surface mesh of wall faces (D124 W1, W2) | **not on main**; was moved to 7t2 | Needed by the Eulerian film only |

## 3. Missing stages and decisions

Numbers are proposals; the next free D at the time of writing is D240.

| ID | What | Needed for | Size (agent time, not sized where marked) |
|---|---|---|---|
| **SP0** | L12 v1 specification page `docs/requirements/step28-sprays.md` with decisions (D240 onward): parcel data layout (D122 (4)), forces, dispersion, injectors, coupling kernel, output, restart, units, verification list; amends L12 | everything | docs only, one PR, 45 to 90 min plus the maintainer's answers to the questions in §7 |
| **SP1** | Particle layer, one rank, one-way: new crate `foehncraft-particles` (SoA, `[[injector]]` keys, parcel ids, tracker driver on `Locator`, sub-stepping by exact exponential integration of the relaxation equation so the step is not limited by the response time, stick and rebound walls, output block, restart block, CPU only) | S2 | 2 to 3 h of code, 1 CI cycle |
| **SP2** | Two-way momentum, dispersion and the first ledger pair: parcel sources deposited per cell, sorted by (cell gid, parcel gid), no atomics (D117 (2)); implicit (S_p) part for stability; stochastic dispersion with the eddy lifetime model (Gosman and Ioannides 1983; Gosman and Ioannides use the k-eps field), optional smoothing through the scalar path (Sun and Xiao 2015); balance line "gas momentum given to parcels = momentum in parcels + wall impulse" | S3, S4 | 2 to 3 h, 1 or 2 CI cycles |
| **SP3** | Parcel migration between ranks: sort parcels by owner after each tracking sub-step, send in a deterministic order, ids stay (D93), load weights for METIS (D122 (1)) | multi-rank sprays | not sized yet: depends on how W3 hands out processor faces; plan after W3 is merged. Do **not** start before W3 |
| **SP4** | Thermodynamic slice of L2: one vapour species transported as a mass fraction on top of 7c (binary mixture: air plus one vapour), mixture density, c_p, mu, k, diffusivity by simple rules (mass-weighted and Wilke 1950 / Fuller et al. 1966 type), ledger quantity "vapour mass"; liquid property set for water, ethanol, n-heptane, n-dodecane as tables from the open literature (NIST webbook values are public domain; check D201 licence entry per source) | SP5 | 3 to 4 h, 2 CI cycles; liquid data entry work separate; full L2 (reaction mechanisms, stiff chemistry) stays out and stays "reacting flows" later |
| **SP5** | Heat and mass transfer, evaporation: Ranz-Marshall (1952) Nusselt and Sherwood numbers; Spalding mass transfer number; Abramzon-Sirignano (1989) film correction; uniform droplet temperature first; sources of vapour, enthalpy and momentum into the gas; ledger pairs for mass and enthalpy; entropy term (D90) | evaporating cases | 3 to 4 h, 2 CI cycles |
| **SP6** | Breakup and collision: TAB (Taylor 1963, O'Rourke and Amsden 1987) or Reitz-Diwakar (1987) bag and stripping first, then KH-RT (Reitz 1987, Beale and Reitz 1999); collision by O'Rourke (1981) sampling or Nordin (2001) | diesel-type and atomiser cases | 3 to 5 h; start only on a decision on which atomisation regime comes first |
| **SP7** | Wall film F1 (D123, D124): surface mesh W1 and W2 (moved to 7t2), Eulerian film with Bai-Gosman regimes | spray impingement | not sized yet; needs 7t2, ray casting (D108 part 2) is optional |
| **SP8** | Rotating zones: parcel velocity in the rotating frame, crossing of sliding interfaces, MRF | mixer, rotary atomiser | not sized yet; leave out of v1 |
| **SP9** | GPU kernels for the parcel step (D47, D60) | speed | not sized yet; CPU first, kernels written as pure functions so the lift is mechanical |

SP2 completes "first useful spray module". SP4 and SP5 complete "evaporating
spray". SP3 can run in parallel with SP4 to SP6 once W3 is merged.

## 4. Order and the smallest first PR

Order: SP0 (spec with answers), SP1, SP2, then SP4, SP5, SP6 in parallel with
SP3 after W3 to W5, then SP7. SP1 can start the moment SP0 is agreed; there is
no code dependency on the MPI work, only the rule that all kernels are pure
functions of (parcel, gas state), so they run unchanged on each rank later.

**Smallest first PR (SP1 slice): "one-way inert parcels in a prescribed flow".**
Parcel layer, one `[[injector]]` (point source, cone, fixed mass flow, droplet
sizes from a Rosin-Rammler or constant diameter), Schiller-Naumann (1935) drag
and gravity, stick and rebound walls, no feedback on the gas, VTKHDF parcel
block, one gallery example, one case file, one recipe, `foehn run` support.
Validation (all analytic or open literature, no code-to-code numbers):

| Check | Reference |
|---|---|
| Terminal velocity of a settling droplet in still air (Stokes regime, then Schiller-Naumann at Re up to 800) | Stokes 1851; Schiller and Naumann 1935; the equation is our own ODE, solved to 1e-10 |
| Velocity relaxation of a parcel released in a uniform stream: v(t) = u(1 - exp(-t/tau_p)) | exact; also proves the exponential sub-stepping |
| Parcel in solid-body rotation or Taylor-Green flow: a fixed streamline-crossing test against an RK4 reference at 1e-3 of the radius | analytic velocity field |
| Point-location contract: every parcel stays in the cell it is assigned to, tetrahedral, polyhedral and hexahedral meshes, after a renumbering | existing `locate.rs` tests, repeated through the tracker |
| Bitwise reproducibility with 1 to 4 threads, parcel ids from (injector, step, index) | D93 |

SP2 then adds the two-way test that decides acceptance: **a closed box with a
falling droplet cloud conserves momentum between gas, parcels and wall to
round-off**, and the settling-cloud drag on a fluid column (Maxey and Riley
1983; Crowe, Schwarzkopf, Sommerfeld, Tsuji 2011, textbook) matches the
buoyant-weight balance of the mixture.

Validation ladder for later stages (all open literature, each checked against
D201 before it goes in a table):

- Particle dispersion in grid turbulence: Snyder and Lumley (1971), Wells and
  Stock (1983).
- d-squared law of an evaporating droplet (Godsave 1953; Spalding 1953), with
  and without convection (Ranz-Marshall; Abramzon-Sirignano 1989).
- Penetration of a non-evaporating spray: Hiroyasu and Arai (1990) correlation,
  Naber and Siebers (1996); Sandia ECN data are public, licence entry to be
  checked first.
- Evaporating droplets in a closed box: mass and energy of gas and parcels sum
  to zero, entropy production nonnegative (already in the L12 text).
- Spray impingement on a flat plate: Katsura et al. (1989), through the wall
  film (SP7).

## 5. Compressibility, reacting flows and species

Compressibility: **not a prerequisite**. 7c (low-Mach variable density with
p0(t) and closed volumes) covers evaporation into a closed or open box,
spray cooling, SCR urea mixers, scrubbers. Fully compressible flow with
shocks or acoustics (engine intake near choking, supersonic injection) is not
planned in the stage list and not part of L12 v1.

Reacting flows: **out of scope** for L12 v1. SP4 adds one passive-like vapour
species transported on the 7c density, with mixture properties; it is a
slice of L2, not L2. Liquid fuel into a reacting flow is L12 plus L2 plus L28
and stays later.

Energy: the gas energy equation exists; the spray adds the sensible heat
exchange per parcel at SP2 (convective heat transfer without phase change, by
Ranz-Marshall) and the latent heat at SP5. That means SP2 already needs the
enthalpy ledger pair, not only SP5.

## 6. Risks

- **Tracking on bad cells.** The walk is robust on convex polyhedra; at
  non-convex cells it falls back to a bin search. The tracker must count
  fallbacks as a health number (D69) and lose no parcel silently (parcel
  loss is an error code).
- **Time-step coupling.** A parcel step shorter than the flow step is free;
  longer is not allowed. With the exponential integrator the step is limited
  only by the cell crossing, not by the response time. The cost per step with
  parcels crossing many cells must be reported in the phase timers.
- **MPI halo.** A parcel can cross more than one halo layer in a step at large
  CFL. The tracker must stop at the rank face and hand the remainder to the
  neighbour (SP3), never jump across more than the halo.
- **Reproducibility.** Parcel sources are the one place where floating-point
  order matters for the flow: sums are by (cell gid, parcel gid); exact sums
  of W2 are not enough because the partial sums are in the deposit.
- **Python and agents.** The spray is a case-file section, a registered tool
  (D196) and an MCP-visible key set from SP1, so the agent-first rule is kept.
- **Clean-room.** Every correlation above is from the open literature (listed
  in §8); no constants taken from the documentation of commercial codes; model
  constants as case keys with the literature value as the default.
- **Rotating zones.** Parcels in MRF or sliding zones need frame handling;
  leave them out of v1 and refuse such a case with a diagnostic code.

## 7. Questions for the maintainer (SP0)

1. **First application.** Which of these decides the validation list and the
   gallery case: spray cooling (sensible heat), urea-water for SCR (film
   needed), diesel-type (breakup and evaporation), painting, scrubber, dust?
   My recommendation: water-spray cooling a hot gas stream, because it needs
   SP1 to SP5 only and has the closed-box ledger test.
2. **Droplet size input.** Constant diameter, Rosin-Rammler, log-normal, or a
   table? Recommendation: all three as keys, constant first.
3. **Dispersion.** Eddy-lifetime model (Gosman and Ioannides) first, or the
   continuous random walk? Recommendation: eddy lifetime.
4. **Ledger.** Do SP2 on the fixed structs (one new struct for the exchange
   pair) or wait for the ledger as data (D119, 7.0a)? Recommendation: do not
   wait; the pair has a fixed form and moves to rows later.
5. **GPU.** Keep CPU only until SP9? Recommendation: yes.
6. **Parcel data in the case file.** Injectors in the case file or in an
   included file (D114)? Recommendation: `[[injector]]` in the case file, the
   table of sizes in an include.
7. **Multi-rank.** SP3 after W3 to W5 are merged, as proposed above, or an
   earlier serial-only release of sprays with a refusal on more than one rank?
   Recommendation: serial release first (SP0 to SP2, SP4, SP5), SP3 later.

## 8. References (open literature)

- Abramzon, B., Sirignano, W. A. (1989), Droplet vaporization model for spray
  combustion calculations, Int. J. Heat Mass Transfer 32.
- Bai, C., Gosman, A. D. (1995), SAE 950283 (see the wall film assessment).
- Beale, J. C., Reitz, R. D. (1999), Modeling spray atomization with the
  Kelvin-Helmholtz/Rayleigh-Taylor hybrid model, Atomization and Sprays 9.
- Crowe, C. T., Schwarzkopf, J. D., Sommerfeld, M., Tsuji, Y. (2011),
  Multiphase Flows with Droplets and Particles, 2nd ed., CRC Press.
- Godsave, G. A. E. (1953), Studies of the combustion of drops in a fuel
  spray, Fourth Symposium (International) on Combustion.
- Gosman, A. D., Ioannides, E. (1983), Aspects of computer simulation of
  liquid-fuelled combustors, J. Energy 7.
- Hiroyasu, H., Arai, M. (1990), Structures of fuel sprays in diesel engines,
  SAE 900475.
- Macpherson, G. B., Nordin, N., Weller, H. G. (2009), Particle tracking in
  unstructured, arbitrary polyhedral meshes for use in CFD and molecular
  dynamics, Commun. Numer. Meth. Engng 25.
- Maxey, M. R., Riley, J. J. (1983), Equation of motion for a small rigid
  sphere in a nonuniform flow, Phys. Fluids 26.
- Naber, J. D., Siebers, D. L. (1996), Effects of gas density and vaporization
  on penetration and dispersion of diesel sprays, SAE 960034.
- Nordin, N. (2001), Complex chemistry modeling of diesel spray combustion,
  PhD thesis, Chalmers.
- O'Rourke, P. J. (1981), Collective drop effects on vaporizing liquid sprays,
  PhD thesis, Princeton.
- O'Rourke, P. J., Amsden, A. A. (1987), The TAB method for numerical
  calculation of spray droplet breakup, SAE 872089.
- Ranz, W. E., Marshall, W. R. (1952), Evaporation from drops, Chem. Eng.
  Prog. 48.
- Reitz, R. D. (1987), Modeling atomization processes in high-pressure
  vaporizing sprays, Atomization and Spray Technology 3.
- Reitz, R. D., Diwakar, R. (1987), Structure of high-pressure fuel sprays,
  SAE 870598.
- Schiller, L., Naumann, A. (1935), Uber die grundlegenden Berechnungen bei
  der Schwerkraftaufbereitung, Z. Ver. Deutsch. Ing. 77.
- Snyder, W. H., Lumley, J. L. (1971), Some measurements of particle velocity
  autocorrelation functions in a turbulent flow, J. Fluid Mech. 48.
- Spalding, D. B. (1953), The combustion of liquid fuels, Fourth Symposium
  (International) on Combustion.
- Sun, R., Xiao, H. (2015), Diffusion-based coarse graining in hybrid
  continuum-discrete solvers, Int. J. Multiphase Flow 72.
- Wells, M. R., Stock, D. E. (1983), The effects of crossing trajectories on
  the dispersion of particles in a turbulent flow, J. Fluid Mech. 136.
