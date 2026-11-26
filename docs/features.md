# What works today

The full list of what Foehncraft does on `main` now. The
[README](../README.md#what-works-today) has the short version; the
[Example gallery](gallery.md) shows every feature with pictures and the
commands, and the [Engineering Gallery](applications.md) puts them together
on industrial geometries.

Everything below is in release 1.0.1 and on `main`, checked by the
verification cases of the [case catalogue](cases.md) and shown in the
[gallery](gallery.md). How well each feature is checked (verified,
partial or unverified, with the evidence and the known gaps) is on the
[verification status](verification/status.md) page (the generated
[verification manual](verification/manual.md) was last run on 0.6.0;
its regeneration is in v1.1), and every run
lists its limitations (a partial feature in use, a value outside the
validated range, one mesh only, …) in its report and its JSON result.

- **Meshes.** Arbitrary polyhedral cells in 2D and 3D; the Gmsh v4, SU2, CGNS (HDF5),
  OpenFOAM polyMesh, Fluent `.msh`, VTK (`.vtk`, `.vtu`) and VTKHDF (also
  blocked) readers behind one mesh-format interface (2D SU2 files extruded
  one layer), a VTU writer;
  a built-in test mesher (hexahedra, tetrahedra, prisms, polyhedra; graded
  segments, the cylinder O-grid, annuli for rotors); a quality report with
  automatic remediation of skewed cells (Gauss pressure gradient, limiter,
  bounded Rhie–Chow face velocity); an operator exactness check that flags
  the cells where the velocity gradient or Hessian is not exact for a
  quadratic field (D236); partitioning with METIS and halos;
  VTKHDF output for ParaView, in blocks that MurmelCraft reads in parts,
  written block by block.
- **Numerics.** Cell-centred finite volumes with weighted least-squares
  gradients, Venkatakrishnan(-Wang) limiters, for the scalars and the
  temperature also the vertex-based MLP limiter with Nishikawa's R3
  (bounded on tetrahedra and polyhedra, D237), linear upwind and bounded
  central convection, non-orthogonal correction, the time-step-independent
  Rhie–Chow interpolation; BDF2 and Crank–Nicolson; reproducible
  (gather-based) assembly in double precision; exact dot products, norms
  and sums on the CPU, bitwise the same for any thread count, cell order
  and number of ranks (MPI step W2); the steady scalar solver, multigrid
  included, and the flow solver (SIMPLEC, PISO and the coupled solver,
  k-ω SST and k-ε, passive scalars and the energy equation) run on several
  ranks with bitwise the result of one (MPI steps W3 to W5); a restart
  gathered from the ranks is byte for byte the one-rank file and continues
  on any rank count (W6); `foehn flow case.toml --ranks N` (or `mpirun -n N
  foehn flow`) solves a run on ranks, steady or transient, laminar or with
  k-ω SST and k-ε, with passive scalars, the energy equation, a rotating
  frame, the rotating zones of `[[mrf]]` and the sliding mesh (conformal
  stepping and imprinting), also in closed vessels, and writes every output
  file byte for byte as one rank does (W7).
- **Pressure–velocity coupling.** SIMPLEC and SIMPLE for steady flow, PISO
  with outer correctors for transient flow (SIMPLEC-type correctors with
  outer iterations, stable on thin wedge cells at an axis), a time step from the
  integration error, and a coupled pressure–velocity algorithm with
  pseudo-time steps and multigrid sweeps (the backward-facing step 6 times
  faster than SIMPLEC).
- **Linear solvers.** Aggregation multigrid with CG for the pressure (mixed
  precision, hypre as the reference), reduction multigrid for momentum and
  turbulence, BiCGStab.
- **Turbulence (RANS).** k-ω SST (2003), standard and realizable k-ε,
  k-ε-ζ-f (Hanjalić, Popovac & Hadžiabdić 2004), k-ε-ζ-f₀ (Popovac 2023)
  and the elliptic-blending k-ε; wall treatments `low-re`, `automatic`,
  `scalable`, `standard`, `compound` and `compound-2023`, rough walls; the
  exact wall distance; a sink of k for the ζ-f models that is ε at the
  wall, so their results converge with the mesh; the curvature and
  rotation correction of SST
  (Smirnov & Menter 2009), a whole domain in a turning frame, and the
  quadratic constitutive relation QCR2000/QCR2013 for the secondary flow
  of non-circular ducts; both converge on tetrahedra and perturbed cells
  through a quadratic least-squares velocity fit (gradient and Hessian)
  on the cells the operator check flags (D236).
  Second order verified with manufactured
  solutions of every model on hex, polyhedral and tetrahedral cells
  (the QCR stress is first order on tetrahedra); validated on the flat plate, the channel against DNS, the backward-facing step and
  the impinging jet.
- **Heat and mass transfer.** Passive scalars (temperature, concentration)
  carried by laminar and turbulent flow, with the turbulent Prandtl number
  (constant or Kays–Crawford) and thermal wall functions (Kader,
  Jayatilleke): the heated pipe, plate and channel against correlations and
  DNS, the Nusselt number of the impinging jet against the experiment.
  The energy equation of the fluid (stage 7): c_p(T), k(T) and μ(T) from
  the material database, viscous heating, temperature, heat-flux, heat-flow
  and convection walls, and a heat balance in the report; second order on
  manufactured and exact solutions. Conjugate heat transfer: solid regions
  of a mesh join the temperature equation, the flow runs on the fluid
  cells and the report gives the heat flow through every fluid–solid
  interface (a composite wall exact to round-off). With a turbulence model
  the thermal wall function sits in series with the solid's half cell at
  every interface, and the turbulent channel over a conducting wall gives
  the same heat flow as the independent wall-function path to 10⁻⁴. In
  time, the fluid and the solids heat up in one matrix with the time terms
  of both: two bodies brought into contact keep the exact contact
  temperature and follow the erf profiles at second order in space and
  time, and a heated wall switched on closes its heat balance (input,
  outflow, heat stored per region, written every step) to 10⁻¹¹ in every
  step. Contact resistances and thin walls sit in the faces between two
  regions, solids may conduct orthotropically (axes per region), a region
  can carry its own heat source, and with one-sided gradients at the
  interfaces a piecewise-linear temperature is exact on a skewed mesh (a
  rotated laminate, a fluid and a block behind a contact resistance, to
  3·10⁻¹¹ K); a chip on a thermal pad on a circuit board with copper
  planes is the gallery example. Thermal acceleration of the solids
  (`[energy.acceleration]`): ρ c_p of the solids is scaled down until they
  settle (the factor measured from each solid's time constant, or given),
  then the true ρ c_p comes back, at once or ramped over a few periods in a
  periodic flow; the heated wall reaches its steady state five times sooner,
  and the heat the true solid would have stored is booked in the balance.
  The cycle-averaged loop (`[energy.cycle]`) couples a periodic flow with
  slow solids: the flow runs a few cycles with the solids held, the solids
  are solved steady under the cycle-mean heat transfer (answering with the
  response of the cycle-mean flow), and the new wall temperatures go back,
  until they settle; a pulsating channel over a heated wall converges in 4
  updates to within 0.4 % of the fully coupled periodic state. With
  `skin = "effusivity"` the held wall swings like a thin skin of its
  effusivity, and `resolved_cycles` ends the loop with a few fully coupled
  cycles: a thin-skin wall ends at 1250 s within 0.016 K and 0.1 % on the
  swing of the periodic state the fully coupled run reaches after about
  5000 s; `solid = "transient"` gives the warm-up history instead, marched
  under the cycle-mean flow with the clock jumped by whole periods (within
  0.6 % of the rise, 780 steps against 6000). The frozen-flow march
  (`[energy.frozen_flow]`) holds a settled flow and marches only the
  temperatures with long steps, refreshing a flow that follows the
  temperature when the fluid has moved: a heated wall under a fluid whose
  viscosity halves over 17 K warms up in 557 steps instead of 6000, within
  0.06 K of the fully coupled run. Conjugate runs, the loop and the march
  restart bit-identically. Several fluids in one case (D194): any number
  of fluid regions, each with its own material (liquids and gases) and
  the temperature at which it is evaluated, beside any number of solids;
  fluids that touch only through solids are separate parts of the flow,
  each with its own pressure level; a water / oil counterflow exchanger
  across a steel plate meets ε–NTU to 0.07 %, and a shell and tube
  exchanger runs hot air, water and steel in one case on 0.48 M
  tetrahedra (application 203, a demonstration on one mesh).
- **Gravity and buoyancy.** `[gravity]` and the Boussinesq buoyancy
  (`[buoyancy] model = "boussinesq"`) in SIMPLEC and the coupled algorithm.
  The solver works with the reduced pressure p′ = p − ρ_ref g·(x − x_ref)
  and writes the static pressure beside it. Values of pressure boundaries
  are p′ by default, or the static pressure with
  `pressure_kind = "static"`. The buoyancy enters the faces like the
  pressure gradient (balanced force), so a stratified box at rest on a
  perturbed hex, tet, poly or mixed mesh stays at rest to 10⁻¹⁵ m/s; the
  buoyant slot meets its exact profile at second order, and a warm strip
  in a closed air enclosure drives a laminar plume (gallery 186). Laminar
  natural convection benchmarks (gallery 272): the heated square cavity
  at Ra 10³ to 10⁶ (mean Nu with GCI ≤ 0.12 %), the vertical plate within
  0.4 % of the Ostrach similarity solution, and mixed convection with
  flow reversal at second order against its exact profile.
- **Variable density.** The low-Mach density ρ(T, p₀) of an ideal gas or a
  liquid in every cell and face (`[energy] density = "variable"`, or
  `[buoyancy] model = "variable-density"` with the body force (ρ − ρ_ref)
  g), with SIMPLEC and the coupled algorithm, steady and transient,
  turbulent and in conjugate heat transfer. p₀ is the operating pressure
  of an open domain, or follows from the fixed mass of a closed volume,
  with the pressure work that makes a closed gas heat with c_v. A heated
  duct meets its exact balances, a closed box heated uniformly its exact
  p₀(t), a manufactured solution with ρ = 1/T converges at second order,
  and air between walls at 960 K and 240 K has a heat flow with GCI
  0.02 % (gallery 195).
- **Lagrangian sprays (L12, SP1).** Droplets as parcels in a transient flow run (`[spray]` with an included injector file): exact exponential integration of drag (Schiller and Naumann) and gravity, tracking on the mesh locator, stick, rebound, escape or mirror at boundaries, a mass book that stops the run when it does not close, parcel files and restart, results bitwise equal on 1 to 4 threads, one source for CPU and GPU kernels. One process ([gallery 242](gallery/sprays.md#spray-duct)).
- **Datums and spray statistics (SP1b).** Rosin-Rammler and log-normal droplet sizes; `[[datum]]` planes (cut through the cells or moved to the faces, discs, rectangles), STL surfaces, boundaries, zone interfaces, boxes, cylinders and zones; the liquid through each datum bitwise equal on any thread count, the gas mass, volume and momentum flow and area averages; D10 to D43, Dv10/50/90, histograms, a Rosin-Rammler fit and the penetration, in `datums.csv`, `spray_hist_<datum>.csv` and `result.json`, also as the tools `post.datums`, `post.droplets` and `post.penetration`; `series.json` lists the files of a run for a viewer ([gallery 265](gallery/sprays.md#spray-datums)).
- **Two-way spray coupling, momentum (SP2a).** `[spray.coupling] mode = "two-way"`: the drag on the droplets is an implicit source of the gas momentum (value plus slope), parcel passes inside the time step until the coupling lag is small (at most 3), a closing step that gives the parcels what the gas received, so the ledger pair cancels to the bit; bitwise equal for any parcel order and thread count, stable and without overshoot at steps of 100 droplet response times. ([gallery 280](gallery/sprays.md#spray-two-way)).
- **Two-way spray coupling, sensible heat (SP2b).** With `[energy]`, every parcel carries a temperature that relaxes exactly towards the gas (Ranz-Marshall), and the heat is an implicit source of the gas energy with the same passes and closing step; the enthalpy pair cancels to the bit. Water droplets cool a 300 °C air stream by about 50 K in the gallery duct ([gallery 283](gallery/sprays.md#spray-cooling)).
- **Turbulent dispersion of the sprays (SP2c).** `[spray] dispersion = "eddy-lifetime"` (Gosman and Ioannides): Gaussian eddies of variance 2k/3 with lifetime and crossing time, keyed per parcel and eddy so runs stay bitwise reproducible on any thread count; tracers spread as the exact random walk to 1.3 % ([gallery 285](gallery/sprays.md#spray-dispersion)).
- **Steady spray mode (SP2s).** In a steady case `[spray.steady]` follows droplet streams through the frozen gas and hands it relaxed momentum and heat sources every few outer iterations; the run is converged only when the gas and the source residual are, the ledger pairs cancel to the bit, and the streams are keyed so runs are bitwise equal on any thread count. The spray-cooled duct converges in 14 s and comes within 3 % of the transient run's heat ([gallery 292](gallery/sprays.md#spray-steady)).
- **Evaporating sprays (SP4, SP5).** `[spray.evaporation]` (Ranz-Marshall or Abramzon-Sirignano) lets water droplets evaporate into a vapour species of `[species]`: d² law and wet-bulb relaxation inside each sub-step, vapour, mass and enthalpy as implicit sources of the gas with the same passes and closing step, the vapour ledger pair cancelling to the bit, and a counted saturation limiter as a safety net. A closed box of air and water ends at the exact adiabatic saturation; 4000 droplets in one nearly saturated cell never oversaturate it and give the same bits on 1 and 4 threads. In the steady mode the duct with more water than the air can take saturates, but its source residual stalls at about 0.1 ([gallery 306](gallery/sprays.md#spray-evaporation)).
- **Spray breakup (SP6).** `[spray.breakup]` breaks fuel jets into drops: blobs from the nozzle hole (`hole_diameter`, `area_contraction`, a measured `rate` table), stripped into child parcels by the nozzle turbulence, the cavitation or the Kelvin-Helmholtz waves, whichever is fastest (KH-ACT, Som and Aggarwal 2010; `kh-rt` without the nozzle terms), Rayleigh-Taylor breakup beyond the breakup length, and the TAB distortion that raises the drag of a drop towards a disc's. The liquid is kept to round-off at every split and the results are the same bits on any thread count. A 60 MPa diesel-like spray reaches the penetration of the Hiroyasu-Arai correlation within 1 % at 1 ms ([gallery 318](gallery/sprays.md#spray-breakup)).
- **High-pressure droplet surface (Spray A, SA1).** `[spray.evaporation] equilibrium = "peng-robinson"` replaces Raoult's law at the droplet surface by the Peng-Robinson phase equilibrium of the fuel with the ambient gases, with the latent heat from the same equation, zero at the critical mixing state (counted parcels). Built-in n-dodecane: vapour pressure within 1 % of measurements up to 655 K, measured N₂ bubble points within 2.7 % on average; at 6 MPa the surface holds 1.5 to 1.9 times Raoult's fuel vapour. The Spray A flow case follows with the breakup model ([gallery 314](gallery/sprays.md#spray-a-surface)).
- **ECN liquid length and vapour penetration (Spray A, SA3).** `[spray.ecn]` writes per injector and step the liquid length of the ECN modelling standard (where the projected liquid volume across the spray falls to 0.2·10⁻³ and 2·10⁻³ mm³/mm²) and the vapour penetration (0.1 % fuel), the same bits on any thread count. In a chamber of still hot air a water spray's liquid length settles after 100 ms while its vapour keeps going ([gallery 319](gallery/sprays.md#spray-ecn)).
- **Adjoint, first stage (L17 A0).** The steady residual R(q) of SIMPLEC and the coupled solver as one callable (`Flow::steady_residual`, steady laminar flow), and the scatter-free reverse sweep of the momentum face kernel from local dual-number Jacobians, checked by a dot-product test and finite differences; the adjoint solve itself follows in A1 ([gallery 271](gallery/verification.md#adjoint-a0)).
- **Porous regions (D270, PM1).** `kind = "porous"` with a Darcy–Forchheimer resistance in the superficial velocity: one value or three principal values along the axes of the region, or the `honeycomb` shorthand for square ceramic channels (Shah & London). A duct filled with the matrix gives the exact pressure drop, a block in a duct reaches it with second order. Axes along the mesh axes, laminar flow, one process, CPU only; with `[energy]` and a `solid`, gas and matrix carry temperatures of their own, coupled by h·a in one energy matrix (D271, PM2; Schumann's bed exact to second order, a honeycomb store charged) ([gallery 273](gallery/porous.md#porous-flow), [275](gallery/porous.md#porous-two-temperatures)).
- **Immersed bodies by penalisation (L13 stage 1).** Cylinders, spheres and boxes put into any mesh without meshing them (`[immersed]`): solid fraction per cell, an implicit Brinkman sink inside the momentum rows (CPU and GPU, SIMPLEC and coupled), the force on every body. Laminar and isothermal only, the wall smeared over a cell: the Re 20 cylinder on a plain box mesh gets drag and Δp within 1 % of the benchmark, lift 35 % high ([gallery 270](gallery/immersed.md#penalised-cylinder)). A density-field topology search (`[immersed.topology]`) designs a duct bend in a few seconds and writes the result as a watertight STL; the penalised equations converge at second order, and a sweep sets the permeability at about 0.01 h² ([gallery 291](gallery/immersed.md#topology-bend)).
- **Immersed bodies as cut cells (L13 stage 2a, 2b).** `[immersed] method = "cut"` cuts the same bodies out of any polyhedral mesh: watertight clipped faces, wall faces fanned around a point on the surface, split and merged small cells, the body surface as an ordinary wall boundary. Volume and area converge at second order on hexahedra, tetrahedra and polyhedra (benchmarked against marching cubes and the dual-contouring surface point); the Re 20 cylinder on a plain box mesh gets drag in the benchmark range, lift 2 % and Δp 1 % below it ([gallery 313](gallery/immersed.md#cut-cylinder)). STL bodies, wall functions, CHT and near-body refinement are next.
- **Deforming meshes and the piston (stage MS0, L22).** `[mesh_motion]` moves the nodes of a transient case by an expression of place and time: swept volumes per face that close the space conservation law to round-off on any cell type, the flux relative to the moving faces, the volume of each time level, walls that move with their nodes. A uniform flow on a swinging mesh stays uniform to 10⁻¹⁵, and a piston compressing a closed gas volume by four follows p₀V^γ with an error that falls as Δt² ([gallery 296](gallery/moving-solids.md#deforming-mesh)). Prescribed motion, one process, CPU; mesh-motion solver and valves moved by the fluid follow (MS1, MS2).
- **Species of gas mixtures (L2, steam furnace stage SF1).** `[species]` with the 53 species of GRI-Mech 3.0: N − 1 mass fractions carried by the flow, the closure species as 1 − Σ, Schmidt or Lewis diffusion, compositions in mole or mass fractions; the density follows the mixture's molar mass on the low-Mach variable density, and c_p, λ, μ are those of the mixture (NASA-7, kinetic theory). Species and element accounts in the ledger. Binary diffusion meets its exact decay, a methane jet in air closes its element accounts to the solver tolerance (gallery [278](gallery/species.md#species-transport)). The species react in the cells with SF2b (below).
- **Gas-phase chemistry (L2, D120, steam furnace stage SF2a).** Reaction mechanisms read from the YAML format of D103 (built in: GRI-Mech 3.0 with 325 reactions, the global methane schemes of Westbrook & Dryer and Jones & Lindstedt): Arrhenius, three-body and Lindemann/Troe fall-off rates, reverse rates from the equilibrium constants, an analytic Jacobian, and the stiff Rosenbrock integrator RODAS3 as a cell-local kernel for CPU and GPU. Chemical equilibrium by Gibbs minimisation; `foehn chemistry reactor` and `foehn chemistry equilibrium`. Methane–air ignites at 3.4 ms from 1400 K and ends 4·10⁻⁵ K from its equilibrium; the adiabatic flame temperature is 2224.9 K (gallery [289](gallery/species.md#gas-chemistry)).
- **Reacting flows in the cells (L2, D120, steam furnace stage SF2b).** `[chemistry]` integrates every cell of a transient laminar flow over each time step with the stiff kernel and feeds the species changes and the heat release back into the flow on the variable density; the same cell body runs on the CPU threads and as a GPU entry (`[run] backend = "cuda"`). A laminar methane–air flame with GRI-Mech 3.0 burns at 0.348 m/s against the measured 0.360 ± 0.017 m/s (gallery [298](gallery/species.md#reacting-flow-cells)). Turbulent and steady reacting flows follow with SF3.
- **Water and steam (D102, steam furnace stage SF7a).** The materials `water-if97` and `steam` take their states from IAPWS-IF97 (regions 1 to 5 and the saturation line) with the IAPWS viscosity, conductivity and surface tension at the operating pressure, for the energy equation and the variable density; `foehn steam` gives any state from (p, T), (p, h), (p, s) or a quality, and the same kernel runs as a GPU entry. The IF97 verification tables are met to the last printed digit, and a heated steam duct meets its exact energy balance (gallery [277](gallery/properties.md#steam-if97)).
- **Internal flows (D45).** Turbulent pipe flow against DNS, the
  asymmetric plane diffuser of Buice & Eaton, the U-bend of Monson et al.
  and a laminar T-junction mixing a passive scalar on polyhedra, each with
  a grid convergence study (GCI) on three meshes and, for the diffuser and
  the U-bend, also with sharp edges (gallery
  [193](gallery/turbulence.md#pipe-turbulent), [196](gallery/flow.md#tjunction-mixing),
  [199](gallery/turbulence.md#diffuser), [204](gallery/turbulence.md#ubend)).
- **Thermal radiation (R1a).** Monte Carlo transport of energy bundles
  through grey absorbing-emitting gas on a prescribed temperature field,
  with grey diffuse walls, symmetry planes and openings: emission from
  scrambled Sobol points, the exact path length in every cell, integer
  tallies with an energy account that closes exactly, bitwise the same on
  any number of threads and ranks; the bundle transport also runs on
  NVIDIA GPUs. Checked against six exact cases (slab, enclosure, parallel
  plates, cube of grey gas) and shown on a box furnace with cooled tubes
  ([gallery 259](gallery/radiation.md#furnace)), where chosen bundles
  of the call are traced and drawn over the net flux into the walls
  ([gallery 281](gallery/radiation.md#furnace-paths)). The coupling
  to the energy equation follows (R1b).
- **Structures (MS5, parts 1 to 4).** A small linear structural finite element
  solver for reed valves and other parts the flow moves (L53):
  Euler–Bernoulli and Timoshenko beams, MITC4 and MITC3+ shells (plates
  clamped along an edge, pressures), quadratic solids (ten-node
  tetrahedra and twenty-node hexahedra from a second-order Gmsh file or a
  generated block, with the von Mises stress), unilateral contact with
  rigid seats and stops, responses in time (generalised-α,
  mode superposition, Rayleigh damping, Craig–Bampton reduction), thermal
  strain loads from temperatures, given or mapped from the solid region
  of a flow case, supports, springs, point masses,
  clamp springs per length of the clamp line (D290) and their fit to a
  measured natural frequency, static deflection and modes by shift-invert
  Lanczos with a Sturm check that no mode was missed, run as the tools
  `structure.run` and `structure.fitclamp`. Checked against cantilever
  deflections and frequencies, the pinned Timoshenko beam and the
  frequency equation of a beam on a rotational clamp spring
  ([gallery 295](gallery/structure.md#reed-beam)), the shells by the
  MacNeal–Harder patch tests, the Navier plate, Leissa's plate modes and
  the Scordelis–Lo roof ([gallery 297](gallery/structure.md#reed-plate)).
  The dynamics are checked against the damped oscillator and the thermal
  loads against free thermal expansion and bending
  ([gallery 299](gallery/structure.md#reed-dynamics)), the solids by
  the patch test, a cantilever against beam theory and Lamé's thick
  cylinder ([gallery 300](gallery/structure.md#reed-stop)), the
  contact by the propped cantilever and a mass bouncing on a stop
  ([gallery 302](gallery/structure.md#reed-valve-contact)), the
  temperature mapping by linear fields mapped exactly and a free body
  that stays stress-free ([gallery 308](gallery/structure.md#heat-sink-hot)).
  The two-way coupling with the flow (FSI) is the rest of MS6.
- **Rotating machinery.** The moving reference frame (steady and transient
  frozen rotor), the sliding mesh by conformal stepping and by imprinting
  the interface in every step, several rotors at their own speeds,
  turbulence on all of them; torque per wall and rotor, shaft power in the
  energy balance, head and hydraulic efficiency; relative and frame
  velocities, flow angles and rothalpy as output.
- **Boundary conditions.** Velocity, mass-flow, static and total pressure
  inlets and outlets, freely combinable and checked; walls (moving,
  turning, with blowing or suction), slip walls and symmetry planes,
  axisymmetric wedges, periodic pairs driven by a mass flow or a pressure
  gradient, inlet profiles from tables, backflow handling, convective
  outlets. Pressure values are absolute static pressures (case format
  0.7, gauge optional): every fluid part takes its operating pressure
  from its own outlet, and the solver keeps the gauge pressure, so an
  absolute case and its gauge twin run the same bits (D239).
  Every region and fluid part starts on its own: start tables per
  region, `"inlet"` temperatures, and each part's start pressure from its
  own boundaries (a Laplace solve between two values, D238).
- **Cases and properties.** TOML case files with checks and guidance;
  boundary values, sources and viscosities as expressions; a material
  database (air, water, oils, glycols) with sources.
- **Runs and results.** Restart (a transient run may also start from a
  steady result), control files (`.stop`, `.result_N`, `.read_toml_N`),
  settings that change on a schedule (also the momentum scheme, from a
  first-order start to linear upwind), probes, forces,
  derived fields (vorticity, Q, λ₂, wall shear stress, y⁺, and the Reynolds stress, scales and wall values of the turbulence models),
  statistics of transient runs (time mean, RMS, covariance, phase averages,
  a confidence interval that says whether the window is long enough), mass and
  kinetic-energy ledgers (named terms per quantity, ready for further quantities), a divergence report, a log that can be
  circulated, JSON records for scripts and agents; an iteration budget
  shown before the run and a diagnosis of a steady run at 2000 iterations
  (trend of each residual, verdict, advice, optional stop), budgets
  per study (cell-hours, GPU-hours, output volume) beside the estimate, and
  the limitations of every run (LIM001 to LIM006, D199).
- **Tools for scripts and agents.** One registry of typed tools (D196):
  `foehn commands` prints the catalogue of every tool with its input and
  output schemas and action kind, `foehn call <id> --input file.json` runs
  any of them from JSON ([tools.md](tools.md)): system, case check,
  estimate, resolve, schema and upgrade, diagnostic codes, materials,
  expressions, meshes, flow runs (streaming their run record, stoppable),
  results, derived fields and the scalar, Poisson and G1 studies. Every
  command and Python function goes through these tools, and `cargo xtask
  layering` checks that no front end reaches the engine past them. The
  Python package reaches the same tools (`foehncraft.call`,
  `foehncraft.commands`); [gallery 171](gallery/runs.md#tool-registry).
  With `--journal FILE` every call goes into a journal with its typed
  input, a run keeps the call that started it in `run.jsonl`, and
  `foehn run --script` replays a journal or a run record without a window
  (D197).
- **An MCP server for assistants.** `foehn mcp` serves every tool to an
  assistant over the Model Context Protocol (L18): only inside the read and
  write roots it is started with, with an action table per project that
  allows, asks the person or refuses each kind of action (start a run,
  exceed the estimate, write outside the roots, overwrite results), run
  limits required, progress and cancellation, summaries with their
  limitations first and histories paged by `result.page`, and every call
  journaled for `foehn run --script`. A running case can be watched and
  steered through its live control port (L19): `run.status` with the live
  residuals, `run.set` for run-time settings (never limits), `run.stop`.
  [Connect an assistant](agents.md), [gallery 175](gallery/runs.md#mcp),
  [gallery 177](gallery/runs.md#control-port). A recorded study of an
  assistant (four fin pitches of the heat sink under run limits and a study
  budget, with the refused calls) replays without a window:
  [gallery 181](gallery/runs.md#heat-sink-study). A recorded chat sizes a
  spray cooler from one question: the assistant writes the case, runs it
  three times and answers with the water flow, the outlet temperature
  across the duct and the limitations ([gallery 316](gallery/runs.md#spray-cooler-chat)).
- **GPU and Python.** The flow loop on NVIDIA GPUs; its results agree
  with the CPU to rounding (bitwise until the exact CPU sums of W2; the
  GPU keeps block sums until v1.1); the Python package `foehncraft` that builds, runs and
  studies cases; the grid convergence index of every mesh sweep.
