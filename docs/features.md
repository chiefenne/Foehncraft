# What works today

The full list of what Foehncraft does on `main` now. The
[README](../README.md#what-works-today) has the short version; the
[Example gallery](gallery.md) shows every feature with pictures and the
commands, and the [Engineering Gallery](applications.md) puts them together
on industrial geometries.

Everything below is in release 1.2.0 and on `main`, checked by the
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
  k-ω SST and k-ε, with passive scalars, species, the energy equation on a
  constant or variable density, statistics, a rotating frame, the rotating
  zones of `[[mrf]]` and the sliding mesh (conformal
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
- **Spray breakup (SP6).** `[spray.breakup]` breaks fuel jets into drops: blobs from the nozzle hole (`hole_diameter`, `area_contraction`, a measured `rate` table), stripped into child parcels by the nozzle turbulence, the cavitation or the Kelvin-Helmholtz waves, whichever is fastest (KH-ACT, Som and Aggarwal 2010; `kh-rt` without the nozzle terms), Rayleigh-Taylor breakup beyond the breakup length, and the TAB distortion that raises the drag of a drop towards a disc's. The liquid is kept to round-off at every split and the results are the same bits on any thread count. A 60 MPa diesel-like spray reaches the penetration of the Hiroyasu-Arai correlation within 1 % at 1 ms ([gallery 318](gallery/sprays.md#spray-breakup)). On ECN Spray G (G1, a 22.5° sector of the 8-hole injector, `wedge` with `dimension = 3`) the untuned model follows the measured liquid and vapour penetration to 0.4 ms and then runs ahead, 33 % (liquid) and 21 % (vapour) at 1.2 ms, with KH-RT as well, while the gas the plumes draw along the axis matches the PIV ([gallery 323](gallery/sprays.md#spray-g)). Since the KH-ACT children have the Kelvin-Helmholtz size of the paper (D348) and with collisions, G1 is 8 % ahead (liquid) and 1 % short (vapour) at 1.2 ms, 23 % ahead at 0.3 ms, with drops of 14.5 µm at 15 mm against 7 to 10 µm ([gallery 344](gallery/sprays.md#spray-child-size)).
- **The drops damp the gas turbulence (D330).** `[spray.coupling] turbulence = "dispersion-work"` takes the work the eddies of the dispersion do on the drops out of k, with the matching ε (or ω) term of KIVA-II, booked only from the part of the drop's motion that follows its eddy, so the term is never positive and free of the noise of the mean slip. The work over ρV is what k loses to 3 % in a lumped cell. Off by default: on ECN Spray G (measured before collisions) it removes nearly all of the plume turbulence and the spray runs further ahead of the measurement ([gallery 333](gallery/sprays.md#spray-turbulence)).
- **High-pressure droplet surface (Spray A, SA1).** `[spray.evaporation] equilibrium = "peng-robinson"` replaces Raoult's law at the droplet surface by the Peng-Robinson phase equilibrium of the fuel with the ambient gases, with the latent heat from the same equation, zero at the critical mixing state (counted parcels). Built-in n-dodecane: vapour pressure within 1 % of measurements up to 655 K, measured N₂ bubble points within 2.7 % on average; at 6 MPa the surface holds 1.5 to 1.9 times Raoult's fuel vapour. The Spray A flow case is stage SA6 ([gallery 314](gallery/sprays.md#spray-a-surface)).
- **Real-fluid droplet film (Spray A, SA2).** With the Peng-Robinson surface the gas film of each droplet takes the Peng-Robinson density and c_p, the dense-fluid viscosity and conductivity factors of Chung et al. (1988) and the high-pressure diffusivity factor of Riazi and Whitson (1993) (`[spray.evaporation.peng_robinson] film = "real-fluid"`, the default). Nitrogen against NIST: factors within 1.3 % (viscosity) and 4.8 % (conductivity) up to 15 MPa. At Spray A conditions the film moves by 1 to 7 % and a droplet lives 3 % shorter; at 15 MPa by up to 14 % ([gallery 324](gallery/sprays.md#spray-a-film)).
- **The Spray A mixing line (SA1 part 2, D319).** `realfluid::Streams` is the adiabatic mixing line of a fuel and an ambient gas at one pressure, frozen, in Peng-Robinson equilibrium (Michelsen's stability test and an isothermal flash) or ideal (ideal gas, Raoult's law), with the saturation mixture fraction Z_sat where the last liquid evaporates; `siebers_liquid_length` gives the mixing-limited liquid length from it. n-dodecane gets its ideal gas (Scott 1974, within 0.11 % of Lemmon and Huber 2004), also as the vapour species `NC12H26` of the flow. At 6 MPa and 900 K the energy balance holds to 7e-16 and Z_sat is 0.357 (equilibrium) against 0.320 (ideal); Siebers' length with these is 20 to 48 % below the measured Spray A lengths, a recorded gap ([gallery 320](gallery/sprays.md#spray-a-mixing-line)).
- **Liquid property library and film state (D320, D321).** `[spray.liquid] library` takes water, n-heptane, n-dodecane or iso-octane with temperature-dependent vapour pressure, density, specific heat, latent heat, surface tension and viscosity (fits to the NIST WebBook tables; iso-octane to published correlations), with the vapour species NC7H16, IC8H18 and NC12H26. The evaporation takes the gas properties of the droplet film by the 1/3 rule from the species data (default `[spray.evaporation] film = "one-third"`), which halves the error against measured n-heptane droplets: +16 to +59 % against Chauveau et al. instead of +93 to +99 % with far-field properties, a gap shared by this model class ([gallery 321](gallery/sprays.md#spray-liquid-library)).
- **ECN liquid length and vapour penetration (Spray A, SA3).** `[spray.ecn]` writes per injector and step the liquid length of the ECN modelling standard (where the projected liquid volume across the spray falls to 0.2·10⁻³ and 2·10⁻³ mm³/mm²) and the vapour penetration (0.1 % fuel), the same bits on any thread count. On a half, quarter or sector domain the projection adds the mirror images of the parcels in the symmetry planes through the axis (#326). In a chamber of still hot air a water spray's liquid length settles after 100 ms while its vapour keeps going ([gallery 319](gallery/sprays.md#spray-ecn)).
- **Distorted drops (Spray A, SA4).** `[spray.breakup] drag = "corrected-distortion"` gives each drop the drag of a viscous spheroid of its TAB distortion and the heat and mass transfer of its shape, larger surface and internal circulation (the Corrected Distortion model of Nguyen, Dahms, Pickett and Tagliante 2022) in place of the disc-like dynamic drag. In an ECN Spray A run at 900 K and 6 MPa it shortens the liquid length from 8.5 to 7.6 mm, the direction of the paper; both stay short of the measured 11.2 mm, a recorded gap ([gallery 326](gallery/sprays.md#spray-a-distortion)).
- **Drop collisions (D331).** `[spray.collision] model = "ntc"` lets the drops of two parcels in one cell collide: the no-time-counter sampling of Schmidt and Rutland (2000) with the outcomes of O'Rourke (1981), coalescence below the critical impact parameter and grazing separation above it, as in the reference set-ups of the ECN. Mass, momentum and liquid enthalpy are kept to round-off, the collision rate matches the exact pair rate within 3 %, and the drops are the same for any thread count or parcel order. In ECN Spray G (G1) it cuts the liquid lead at 0.8 ms from 19 % to 12 % ([gallery 335](gallery/sprays.md#spray-collisions)).
- **Spray injection from measured droplets (D322) and the Brno spray.** `[[spray.injector]] samples = "table.csv"` starts every parcel or steady stream from one measured droplet of a plane (radius, diameter, axial, radial and optional tangential velocity, optional mass weight), drawn by mass and placed on a random azimuth, and `sector` with `flow_share` models part of a round spray between symmetry planes. The case `spray_brno` runs the air-blast water spray of Rácz et al. from its measured plane at z = 20 mm and compares the planes at 40 and 60 mm per size class with the phase Doppler data ([gallery 325](gallery/sprays.md#spray-brno), recipe 325).
- **Injector devices (D340 to D345).** `[[spray.device]]` gives an injector as it is mounted (tip, axis, zero of the azimuth, clocking) and drilled (`[[device.hole]]`, or rings of equal holes in `[[device.pattern]]`); an `[[injector]]` may name `device` and `hole` instead of a position and direction, `hole = "all"` makes one injector per hole, and old injector files read exactly as before ([gallery 341](gallery/sprays.md#spray-device)). The device carries a surface, a standard gasoline or diesel tip of our own drawing or the user's STL in its own frame, written once per run as `spray/devices.vtu` with the hole axes and injector cones and listed in `series.json`; a picture, not a wall ([gallery 343](gallery/sprays.md#spray-device-surface)).
- **Adjoint, first stage (L17 A0).** The steady residual R(q) of SIMPLEC and the coupled solver as one callable (`Flow::steady_residual`, steady laminar flow), and the scatter-free reverse sweep of the momentum face kernel from local dual-number Jacobians, checked by a dot-product test and finite differences; the adjoint solve itself follows in A1 ([gallery 271](gallery/verification.md#adjoint-a0)).
- **Porous regions (D270, PM1).** `kind = "porous"` with a Darcy–Forchheimer resistance in the superficial velocity: one value or three principal values along the axes of the region, or the `honeycomb` shorthand for square ceramic channels (Shah & London). A duct filled with the matrix gives the exact pressure drop, a block in a duct reaches it with second order. Axes along the mesh axes, laminar flow, one process, CPU only; with `[energy]` and a `solid`, gas and matrix carry temperatures of their own, coupled by h·a in one energy matrix (D271, PM2; Schumann's bed exact to second order, a honeycomb store charged) ([gallery 273](gallery/porous.md#porous-flow), [275](gallery/porous.md#porous-two-temperatures)).
- **Two fluids by the volume of fluid (L3, MP1).** `[vof]` puts water and air, or any two immiscible fluids, into a transient laminar case. α is advected sharply by THINC (or CICSAM), gravity and surface tension act as balanced face forces, and the pressure rows keep the volume. It runs in the same device kernels as the single-phase solver. A water column collapsing in air follows the measured front of Martin & Moyce (1952) to within 5 to 17 %, independent of the mesh. Zalesak's disk keeps its slot (L1 0.048), and the water volume is conserved to round-off. Isothermal, one process; energy and phase change come next (MP2) ([gallery 317](gallery/multiphase.md#vof-dam-break)).
- **Immersed bodies by penalisation (L13 stage 1).** Cylinders, spheres and boxes put into any mesh without meshing them (`[immersed]`): solid fraction per cell, an implicit Brinkman sink inside the momentum rows (CPU and GPU, SIMPLEC and coupled), the force on every body. Laminar and isothermal only, the wall smeared over a cell: the Re 20 cylinder on a plain box mesh gets drag and Δp within 1 % of the benchmark, lift 35 % high ([gallery 270](gallery/immersed.md#penalised-cylinder)). A density-field topology search (`[immersed.topology]`) designs a duct bend in a few seconds and writes the result as a watertight STL; the penalised equations converge at second order, and a sweep sets the permeability at about 0.01 h² ([gallery 291](gallery/immersed.md#topology-bend)).
- **Immersed bodies as cut cells (L13 stage 2a, 2b).** `[immersed] method = "cut"` cuts the same bodies out of any polyhedral mesh: watertight clipped faces, wall faces fanned around a point on the surface, split and merged small cells, the body surface as an ordinary wall boundary. Volume and area converge at second order on hexahedra, tetrahedra and polyhedra (benchmarked against marching cubes and the dual-contouring surface point); the Re 20 cylinder on a plain box mesh gets drag in the benchmark range, lift 1.4 % and Δp 1 % below it ([gallery 313](gallery/immersed.md#cut-cylinder)). CHT and near-body refinement are next.
- **Wall functions on cut bodies (L13 stage 2c).** A turbulent case with wall functions on a cut body: every cut and merged wall cell takes the exact distance of its centroid to the body (D108), the wall function samples the velocity one cell layer out along the wall normal (`[immersed.cut] wall_sample`), the cut wall faces carry the wall shear stress without a non-orthogonal correction, and LIM011 gives the y⁺ range of every body. A flat plate tilted by 2° against the mesh (`shape = "plane"`) gets C_f within 7 % of White's correlation from s = 0.5 on, single faces within about ±15 % (±45 % without the sample) ([gallery 357](gallery/immersed.md#cut-plate)).
- **Conjugate heat transfer through cut cells (L13 stage 2d).** A cut body named as the zone of a solid region (`[regions.<name>] zones = ["<body>"]`) is kept as cells of that region: the cut takes the same crossing points and cut faces for both sides, the fluid and solid parts of every cut cell share each cut face, and the shared faces are the region interface (D160) with the series conductance; small solid pieces are merged into solid neighbours, bodies thinner than about two cells are refused (IMM205, pointing to the thin wall of D189). Concentric cylinders converge at second order for conductivity ratios 9 to 900 with the heat balance to round-off; a heated pin in the Schäfer–Turek channel ([gallery 359](gallery/immersed.md#cut-cht)).
- **Conditions on cut bodies (L13 stage 2e).** Named surface patches of an STL body (`patches`: the solids of an ASCII file or the attribute tags of a binary one) become boundaries `<body>-<patch>` with any condition (inlets, outlets, walls); `omega` turns the wall of a body about its axis with the torque in the boundary summary; `invert` keeps the fluid inside a cylinder, sphere or box. Taylor–Couette flow with both cylinders cut out of a box mesh: velocity second order, torque within 0.23 % on 176 cells across; a cylinder blowing through a slot ([gallery 360](gallery/immersed.md#cut-patches)).
- **Near-body refinement with hanging-node faces (L13 stage 2f).** `[immersed.refine]` halves the hexahedra near the immersed bodies up to `levels` times before the cut (`layers` cells per level, or `distance` in metres), 2:1 balanced across faces, edges and corners; the faces between levels are split into the faces of the finer side with the hanging nodes on every edge, so the solver sees a closed polyhedral mesh (D57). Volume and face areas are kept to round-off, the result is deterministic and serial; planar meshes are refined in their plane, other cell types kept whole. A conduction MMS converges at second order across the hanging faces; the Schäfer–Turek cylinder on a D/10 box refined to D/40 near the body matches the uniform D/40 mesh with a seventh of its cells (IMM206; [gallery 361](gallery/immersed.md#refine)).
- **STL bodies (L13 stage 2).** `shape = "stl"` reads a triangle surface (`file`, `scale`) for the penalisation (exact cut volumes as solid fractions) and for the cut: exact segment–triangle crossings, the vertex of dual contouring where a sharp edge or corner runs through a cell, face points that bend the clipped faces onto the surface (`[immersed.cut] face_point`), a merge that keeps every face valid and the cell volumes exact (faces at a bend stay separate faces, so two cells may share several, which the solver treats as one face) and a per-cell fallback to fewer points where the discretisation would refuse a face. Against the exact moments of the triangles the cut converges at second order on hexahedra, tetrahedra and polyhedra, also in the wall area of a sharp box (first order without the dual-contouring point); Stokes flow through an array of STL spheres reaches the series within 0.9 % on 48³ cells, second order ([gallery 322](gallery/immersed.md#stl-cut)).
- **Deforming meshes and the piston (stage MS0, L22).** `[mesh_motion]` moves the nodes of a transient case by an expression of place and time: swept volumes per face that close the space conservation law to round-off on any cell type, the flux relative to the moving faces, the volume of each time level, walls that move with their nodes. A uniform flow on a swinging mesh stays uniform to 10⁻¹⁵, and a piston compressing a closed gas volume by four follows p₀V^γ with an error that falls as Δt² ([gallery 296](gallery/moving-solids.md#deforming-mesh)). Prescribed motion; on ranks and from restart files bitwise as one rank and as the unbroken run ([gallery 307](gallery/moving-solids.md#deforming-mesh-ranks)); mesh-motion solver and valves moved by the fluid follow (MS1, MS2).
- **Species of gas mixtures (L2, steam furnace stage SF1).** `[species]` with the 53 species of GRI-Mech 3.0: N − 1 mass fractions carried by the flow, the closure species as 1 − Σ, Schmidt or Lewis diffusion, compositions in mole or mass fractions; the density follows the mixture's molar mass on the low-Mach variable density, and c_p, λ, μ are those of the mixture (NASA-7, kinetic theory). Species and element accounts in the ledger. Binary diffusion meets its exact decay, a methane jet in air closes its element accounts to the solver tolerance (gallery [278](gallery/species.md#species-transport)). The species react in the cells with SF2b (below).
- **Gas-phase chemistry (L2, D120, steam furnace stage SF2a).** Reaction mechanisms read from the YAML format of D103 (built in: GRI-Mech 3.0 with 325 reactions, the global methane schemes of Westbrook & Dryer and Jones & Lindstedt): Arrhenius, three-body and Lindemann/Troe fall-off rates, reverse rates from the equilibrium constants, an analytic Jacobian, and the stiff Rosenbrock integrator RODAS3 as a cell-local kernel for CPU and GPU. Chemical equilibrium by Gibbs minimisation; `foehn chemistry reactor` and `foehn chemistry equilibrium`. Methane–air ignites at 3.4 ms from 1400 K and ends 4·10⁻⁵ K from its equilibrium; the adiabatic flame temperature is 2224.9 K (gallery [289](gallery/species.md#gas-chemistry)).
- **Reacting flows in the cells (L2, D120, steam furnace stage SF2b).** `[chemistry]` integrates every cell of a transient laminar flow over each time step with the stiff kernel and feeds the species changes and the heat release back into the flow on the variable density; the same cell body runs on the CPU threads and as a GPU entry (`[run] backend = "cuda"`). A laminar methane–air flame with GRI-Mech 3.0 burns at 0.348 m/s against the measured 0.360 ± 0.017 m/s (gallery [298](gallery/species.md#reacting-flow-cells)).
- **Turbulent flames (L2, D349, steam furnace stage SF3).** `[chemistry] interaction = "edc"` burns the fine structures of every cell over their residence time with the stiff kernel (eddy dissipation concept, Magnussen 2005 constants, version 1981 selectable), `"edm"` is the eddy dissipation model with an Arrhenius cap as the fallback; steady and transient RANS, heat release into the energy equation, failed and guarded cells counted (CHM004). A piloted methane jet flame with the inlets of Sandia Flame D runs steady with both (gallery [351](gallery/species.md#turbulent-flames)); the comparison with the measurements is a gap.
- **Water and steam (D102, steam furnace stage SF7a).** The materials `water-if97` and `steam` take their states from IAPWS-IF97 (regions 1 to 5 and the saturation line) with the IAPWS viscosity, conductivity and surface tension at the operating pressure, for the energy equation and the variable density; `foehn steam` gives any state from (p, T), (p, h), (p, s) or a quality, and the same kernel runs as a GPU entry. The IF97 verification tables are met to the last printed digit, and a heated steam duct meets its exact energy balance (gallery [277](gallery/properties.md#steam-if97)).
- **Water and steam circuits (L3, steam furnace stage SF7b).** One-dimensional tube circuits of a boiler on the IAPWS-IF97 states: tubes, headers, pressure nodes and steam drums, solved together by a Newton method, so that the natural circulation of a drum loop follows from gravity against the losses. Two-phase friction by Friedel or the homogeneous model, void by drift flux (Rouhani–Axelsson), wall temperatures by Gnielinski, Chen or Gungor–Winterton, subcooled boiling from the Saha–Zuber point with the Levy profile, the critical heat flux margin by Bowring or a look-up table, run as the tool `circuits.run`. Checked by the exact split of laminar parallel tubes, the heat and momentum balances of a riser and the closed mass, energy and pressure balances of a drum loop ([gallery 312](gallery/circuits.md#water-steam-circuits)). The coupling to the 3D furnace walls follows below (SF7c).
- **Circuits behind the walls (steam furnace stage SF7c, D350).** `[circuits]` puts a circuit behind the wall boundaries of a steady flow with the energy equation. Each `[[circuits.walls]]` entry maps the faces of one wall to the segments of one tube by their position along the tube axis. The water side enters the 3D energy matrix as a Robin condition per face, U (T_water − T_b), with U the series of `wall_resistance`, `fouling` and the inner heat transfer coefficient of the segment. U and T_water are linear along the tube between the segment midpoints. Every `every` outer iterations the face heats are summed per segment (exact sums) and the circuit is solved again from its last solution. The run converges when the flow and the coupling residual Σ|Q_3D − Q_circuit|/Σ|Q_3D| both meet their tolerances. The ledger books the circuit's thermal account and its exchange pair with the wall conduction. Checked by a slab on a water tube against an independent march of the heat balance (outlet within 1.2·10⁻⁴ K), the same slab on a boiling tube, and a convective evaporator pass whose two walls share one drum loop ([gallery 356](gallery/circuits.md#circuit-walls)). Steady runs on one rank (CPU, GPU); transient circuits are SF18.
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
  derived fields (vorticity, Q, λ₂, wall shear stress, y⁺, the wall heat flux, heat transfer coefficient and Nusselt number of the energy equation, and the Reynolds stress, scales and wall values of the turbulence models),
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
