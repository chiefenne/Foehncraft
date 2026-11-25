# Example gallery

What Foehncraft computes, one picture per case, with the commands that
reproduce it. Command sequences for everyday tasks are in
[Recipes](recipes.md); all case files are listed in the
[Case catalogue](cases.md).

All commands run from the repository root. Case files are in
`examples/cases/`.

The gallery is split by topic; each page has its own table of contents. The numbers are global and stay the same, so "gallery 12" always means the same entry.

| Page | What it covers | Examples |
|---|---|---|
| [Kernel benchmarks and GPU runs](gallery/benchmarks.md) | gate G1 on CPUs and GPUs, on box, polyhedral and tetrahedral meshes, and the flow loop on a GPU | 1, 2, 7, 16 |
| [Meshes](gallery/meshes.md) | the built-in test mesher in 2D and 3D, partitioning, the Gmsh, SU2, OpenFOAM, CGNS, Fluent `.msh` and VTK readers, flow on poor Gmsh meshes with remediation, the operator exactness check, the canonical local mesh and exact reductions over ranks | 3, 4, 5, 6, 26, 127, 129, 134, 138, 189, 200, 219, 230, 231, 234, 244, 262, 279, 290, 301 |
| [Scalar transport and the pressure solver](gallery/scalar-pressure.md) | convection–diffusion order of accuracy, convection schemes and limiters, the pressure multigrid | 8, 9, 10, 11 |
| [Incompressible flow](gallery/flow.md) | steady and transient laminar verification and benchmarks (cavity, Poiseuille, Taylor–Green, cylinder, time-step control), the coupled algorithm and its multigrids, heat transfer in turbulent flow, rotating machinery (moving reference frame, sliding mesh, stirred vessels, several rotors) and the T-junction mixing study | 12, 13, 14, 15, 20, 78, 80, 81, 83, 85, 86, 87, 94, 95, 96, 98, 99, 100, 101, 102, 225, 103, 196 |
| [Passive scalars and the wall distance](gallery/scalars.md) | passive scalars carried by the flow, the bounded BDF2, the exact wall distance, the MLP vertex limiter | 63, 211 |
| [Energy equation and heat transfer](gallery/energy.md) | the energy equation of the fluid with temperature-dependent properties, viscous heating and thermal walls: a heated cavity, conduction with k(T), Couette flow heated by friction, conjugate heat transfer with solid regions, laminar and turbulent (a finned heat sink), transient, contact resistances, thermal acceleration, the cycle-averaged loop and the frozen-flow march, two fluids in one case, start values per region, gravity with Boussinesq buoyancy and the low-Mach variable density | 154, 156, 158, 159, 160, 180, 190, 202, 208, 206, 186, 195, 272, 221, 238 |
| [Thermal radiation](gallery/radiation.md) | Monte Carlo radiation of grey gas and grey walls on a prescribed temperature field: the box furnace with cooled tubes, the cases with exact solutions, scrambled Sobol against pseudo-random emission | 259, 281 |
| [Turbulence models](gallery/turbulence.md) | the RANS turbulence models with their wall treatments (k-ω SST, the k-ε family, k-ε-ζ-f, k-ε-ζ-f₀, elliptic-blending k-ε) on the flat plate, the channel, the backward-facing step and the impinging jet; curvature correction and QCR; manufactured solutions; the internal flows of D45 (pipe, plane diffuser, U-bend) with GCI | 68, 71, 73, 74, 88, 89, 90, 113, 115, 141, 140, 184, 194, 197, 193, 199, 204 |
| [Boundary conditions](gallery/boundary-conditions.md) | slip walls, axisymmetric wedges, inlet profiles from tables, backflow at pressure outlets, convective outlets, periodic modules, absolute pressure values (D239) | 21, 67, 220 |
| [Lagrangian sprays](gallery/sprays.md) | droplets as Lagrangian parcels in the gas flow: injectors, drag, gravity, wall behaviour, the mass book | 242, 265, 280, 283, 285, 292 |
| [Porous regions](gallery/porous.md) | flow through porous matrices with a direction-dependent Darcy–Forchheimer resistance: a ceramic honeycomb, the exact pressure drop of a block | 273, 275 |
| [Structures](gallery/structure.md) | the linear structural solver for reeds and other parts the flow moves: beams and MITC shells, clamp springs and their fit to a measured frequency, static, modal and transient analyses, thermal loads (MS5) | 295, 297, 299 |
| [Immersed bodies](gallery/immersed.md) | bodies in the flow without a body-fitted mesh: Brinkman penalisation of primitive bodies (L13 stage 1) | 270, 291 |
| [Moving meshes and solids](gallery/moving-solids.md) | meshes whose nodes move in time: the deforming mesh with swept volumes, a piston compressing a closed gas volume, free-stream preservation (stage MS0, L22) | 296 |
| [Species and reacting flows](gallery/species.md) | species transport of ideal-gas mixtures on the variable density: mixture properties from kinetic theory, the closure species, the element ledger (steam furnace path, SF1) | 278, 289 |
| [Expressions and materials](gallery/properties.md) | boundary values, body forces and viscosities as expressions, the material database, regions, walls with a mass flux | 23 |
| [Python](gallery/python.md) | cases built, run and studied from Python: expressions with Python operators, callbacks, study tables | 24 |
| [Verification](gallery/verification.md) | the discretisation uncertainty of mesh sweeps (grid convergence index) and the generated verification manual | 25, 170, 271 |
| [Post-processing](gallery/post-processing.md) | derived fields: vorticity, Q, λ₂, wall shear stress, y⁺, from the run or a restart file | 22, 108, 110 |
| [Run control, health and adjoint readiness](gallery/runs.md) | adjoint-ready kernels, runs that a program drives, balances and the divergence report, the log, the table columns, scheduled changes, iteration budgets, the tool registry and journal, the MCP server, the control port and an assistant's study | 17, 18, 19, 27, 28, 29, 166, 171, 172, 175, 177, 181 |

## All examples by number

| # | Example | Page |
|---|---|---|
| <a id="g1-cpu"></a>1 | [Kernel benchmark G1 on the CPU](gallery/benchmarks.md#g1-cpu) | [Kernel benchmarks and GPU runs](gallery/benchmarks.md) |
| <a id="g1-gpu"></a>2 | [Kernel benchmark G1 on a Tesla T4](gallery/benchmarks.md#g1-gpu) | [Kernel benchmarks and GPU runs](gallery/benchmarks.md) |
| <a id="mesh-types-3d"></a>3 | [The four test-mesh types in 3D](gallery/meshes.md#mesh-types-3d) | [Meshes](gallery/meshes.md) |
| <a id="poly-partition"></a>4 | [A polyhedral mesh in four parts](gallery/meshes.md#poly-partition) | [Meshes](gallery/meshes.md) |
| <a id="mesh-types-2d"></a>5 | [The four test-mesh types in 2D](gallery/meshes.md#mesh-types-2d) | [Meshes](gallery/meshes.md) |
| <a id="tjunction"></a>6 | [T-junction from Gmsh](gallery/meshes.md#tjunction) | [Meshes](gallery/meshes.md) |
| <a id="g1-stage1"></a>7 | [G1 on polyhedral and tetrahedral meshes](gallery/benchmarks.md#g1-stage1) | [Kernel benchmarks and GPU runs](gallery/benchmarks.md) |
| <a id="scalar-mms"></a>8 | [Convection–diffusion, manufactured solution](gallery/scalar-pressure.md#scalar-mms) | [Scalar transport and the pressure solver](gallery/scalar-pressure.md) |
| <a id="oblique-step"></a>9 | [Oblique step: five convection schemes](gallery/scalar-pressure.md#oblique-step) | [Scalar transport and the pressure solver](gallery/scalar-pressure.md) |
| <a id="poisson-iterations"></a>10 | [Pressure solver: iterations over the mesh size](gallery/scalar-pressure.md#poisson-iterations) | [Scalar transport and the pressure solver](gallery/scalar-pressure.md) |
| <a id="poisson-amg-mms"></a>11 | [Poisson manufactured solution with the multigrid](gallery/scalar-pressure.md#poisson-amg-mms) | [Scalar transport and the pressure solver](gallery/scalar-pressure.md) |
| <a id="cavity"></a>12 | [Lid-driven cavity against Ghia et al.](gallery/flow.md#cavity) | [Incompressible flow](gallery/flow.md) |
| <a id="poiseuille"></a>13 | [Poiseuille flow: order of accuracy of the flow solver](gallery/flow.md#poiseuille) | [Incompressible flow](gallery/flow.md) |
| <a id="taylor-green"></a>14 | [Taylor–Green vortex: second order in time and space](gallery/flow.md#taylor-green) | [Incompressible flow](gallery/flow.md) |
| <a id="cylinder"></a>15 | [Cylinder in a channel: the Schäfer & Turek benchmark](gallery/flow.md#cylinder) | [Incompressible flow](gallery/flow.md) |
| <a id="gpu-flow"></a>16 | [The flow loop on a GPU](gallery/benchmarks.md#gpu-flow) | [Kernel benchmarks and GPU runs](gallery/benchmarks.md) |
| <a id="adjoint-ready"></a>17 | [Adjoint-ready kernels: derivatives and checkpointing](gallery/runs.md#adjoint-ready) | [Run control, health and adjoint readiness](gallery/runs.md) |
| <a id="run-contract"></a>18 | [A run that a program can drive](gallery/runs.md#run-contract) | [Run control, health and adjoint readiness](gallery/runs.md) |
| <a id="health"></a>19 | [Balances and the divergence report](gallery/runs.md#health) | [Run control, health and adjoint readiness](gallery/runs.md) |
| <a id="time-control"></a>20 | [Time step from the error, output at exact times, potential-flow start](gallery/flow.md#time-control) | [Incompressible flow](gallery/flow.md) |
| <a id="boundary-conditions"></a>21 | [Boundary conditions: slip walls, wedges, tables, backflow, convective outlets](gallery/boundary-conditions.md#boundary-conditions) | [Boundary conditions](gallery/boundary-conditions.md) |
| <a id="derived-fields"></a>22 | [Derived fields: vorticity, Q, λ₂, wall shear stress and y⁺](gallery/post-processing.md#derived-fields) | [Post-processing](gallery/post-processing.md) |
| <a id="expressions-materials"></a>23 | [Expressions and materials: variable viscosity, blowing and suction, water from the database](gallery/properties.md#expressions-materials) | [Expressions and materials](gallery/properties.md) |
| <a id="python-study"></a>24 | [A parametric study from Python: Couette flow with a viscosity gradient](gallery/python.md#python-study) | [Python](gallery/python.md) |
| <a id="gci"></a>25 | [Grid convergence index: Couette flow and the cylinder at Re = 20](gallery/verification.md#gci) | [Verification](gallery/verification.md) |
| <a id="skewed-tets"></a>26 | [Flow on raw Gmsh tetrahedra: the cell-quality remediation](gallery/meshes.md#skewed-tets) | [Meshes](gallery/meshes.md) |
| <a id="log"></a>27 | [A log that can be circulated](gallery/runs.md#log) | [Run control, health and adjoint readiness](gallery/runs.md) |
| <a id="log-columns"></a>28 | [A table with the columns you choose](gallery/runs.md#log-columns) | [Run control, health and adjoint readiness](gallery/runs.md) |
| <a id="schedule"></a>29 | [Settings that change during the run](gallery/runs.md#schedule) | [Run control, health and adjoint readiness](gallery/runs.md) |
| <a id="passive-scalars"></a>63 | [A heated channel and a bounded front](gallery/scalars.md#passive-scalars) | [Passive scalars and the wall distance](gallery/scalars.md) |
| <a id="periodic"></a>67 | [Periodic modules: a channel and a short heated pipe](gallery/boundary-conditions.md#periodic) | [Boundary conditions](gallery/boundary-conditions.md) |
| <a id="sst"></a>68 | [k-ω SST on the TMR flat plate and in a channel](gallery/turbulence.md#sst) | [Turbulence models](gallery/turbulence.md) |
| <a id="k-epsilon"></a>71 | [k-ε family with wall functions: flat plate, channel, Moody chart](gallery/turbulence.md#k-epsilon) | [Turbulence models](gallery/turbulence.md) |
| <a id="zeta-f"></a>73 | [k-ε-ζ-f to the wall: channel and flat plate](gallery/turbulence.md#zeta-f) | [Turbulence models](gallery/turbulence.md) |
| <a id="zeta-f0"></a>74 | [k-ε-ζ-f₀ at any y⁺: channel, flat plate, backward-facing step](gallery/turbulence.md#zeta-f0) | [Turbulence models](gallery/turbulence.md) |
| <a id="coupled"></a>78 | [Coupled algorithm and pseudo-time steps: the step case in 65 iterations](gallery/flow.md#coupled) | [Incompressible flow](gallery/flow.md) |
| <a id="momentum-multigrid"></a>80 | [A multigrid for the momentum sweep of the coupled algorithm](gallery/flow.md#momentum-multigrid) | [Incompressible flow](gallery/flow.md) |
| <a id="turbulence-multigrid"></a>81 | [A multigrid for the turbulence equations: the step settles in 40 and 65 iterations](gallery/flow.md#turbulence-multigrid) | [Incompressible flow](gallery/flow.md) |
| <a id="pressure-hierarchy"></a>83 | [The pressure hierarchy follows its matrix: the fine step 6 times faster than SIMPLEC](gallery/flow.md#pressure-hierarchy) | [Incompressible flow](gallery/flow.md) |
| <a id="thermal-wall-function"></a>85 | [Heat transfer in turbulent flow: the heated pipe and the heated flat plate](gallery/flow.md#thermal-wall-function) | [Incompressible flow](gallery/flow.md) |
| <a id="heated-channel-dns"></a>86 | [The heated channel against DNS: Nu at Re_τ = 395](gallery/flow.md#heated-channel-dns) | [Incompressible flow](gallery/flow.md) |
| <a id="channel-dns-reynolds"></a>87 | [The heated channel at every DNS point: Re_τ 180 to 1020, Pr 0.71 and 0.025](gallery/flow.md#channel-dns-reynolds) | [Incompressible flow](gallery/flow.md) |
| <a id="elliptic-blending"></a>88 | [Elliptic-blending k-ε: channel at three Re_τ, flat plate, backward-facing step](gallery/turbulence.md#elliptic-blending) | [Turbulence models](gallery/turbulence.md) |
| <a id="step-validation"></a>89 | [Backward-facing step against the experiment: six models, three meshes](gallery/turbulence.md#step-validation) | [Turbulence models](gallery/turbulence.md) |
| <a id="impinging-jet"></a>90 | [Impinging jet at H/D = 2: a pipe profile as the inlet, Nusselt number, the stagnation-point anomaly](gallery/turbulence.md#impinging-jet) | [Turbulence models](gallery/turbulence.md) |
| <a id="mrf"></a>94 | [Rotating machinery, steady: the moving reference frame on Taylor–Couette flow](gallery/flow.md#mrf) | [Incompressible flow](gallery/flow.md) |
| <a id="sliding-mesh"></a>95 | [Rotating machinery, transient: the sliding mesh on Taylor–Couette flow](gallery/flow.md#sliding-mesh) | [Incompressible flow](gallery/flow.md) |
| <a id="sliding-mesh-imprint"></a>96 | [The sliding mesh with a pitch ratio that is not one: the interface imprinted in every step](gallery/flow.md#sliding-mesh-imprint) | [Incompressible flow](gallery/flow.md) |
| <a id="stirred-vessel"></a>98 | [A stirred vessel: blades and baffles with the moving reference frame and the sliding mesh](gallery/flow.md#stirred-vessel) | [Incompressible flow](gallery/flow.md) |
| <a id="rotor-fields"></a>99 | [The flow a blade sees: relative and frame velocity, flow angles and rothalpy](gallery/flow.md#rotor-fields) | [Incompressible flow](gallery/flow.md) |
| <a id="turbulent-rotor"></a>100 | [A turbulent rotor: the stirred vessel with water, in the frame and on the sliding mesh](gallery/flow.md#turbulent-rotor) | [Incompressible flow](gallery/flow.md) |
| <a id="torque-power"></a>101 | [Torque, shaft power and the performance of a machine: the rotor of the stirred vessel](gallery/flow.md#torque-power) | [Incompressible flow](gallery/flow.md) |
| <a id="several-rotors"></a>102 | [Several rotors at their own speeds, and the transient frozen rotor](gallery/flow.md#several-rotors) | [Incompressible flow](gallery/flow.md) |
| <a id="revolved-interfaces"></a>103 | [Interfaces of any rotationally symmetric shape: a disc impeller in a pot, a turning duct section between cones](gallery/flow.md#revolved-interfaces) | [Incompressible flow](gallery/flow.md) |
| <a id="turbulence-derived-fields"></a>108 | [Derived fields of the turbulence models: Reynolds stress, scales, wall values](gallery/post-processing.md#turbulence-derived-fields) | [Post-processing](gallery/post-processing.md) |
| <a id="transient-statistics"></a>110 | [Statistics of a transient run: the mean and the fluctuations of a vortex street](gallery/post-processing.md#transient-statistics) | [Post-processing](gallery/post-processing.md) |
| <a id="curvature-correction"></a>113 | [Curvature and rotation correction of SST: channel in spanwise rotation, backward-facing step](gallery/turbulence.md#curvature-correction) | [Turbulence models](gallery/turbulence.md) |
| <a id="qcr"></a>115 | [Quadratic constitutive relation of SST: secondary flow in a square duct](gallery/turbulence.md#qcr) | [Turbulence models](gallery/turbulence.md) |
| <a id="su2-airfoil"></a>127 | [A PyAero airfoil mesh from an SU2 file](gallery/meshes.md#su2-airfoil) | [Meshes](gallery/meshes.md) |
| <a id="foam-cgns-polyhedra"></a>129 | [A polyhedral mesh from OpenFOAM and CGNS](gallery/meshes.md#foam-cgns-polyhedra) | [Meshes](gallery/meshes.md) |
| <a id="five-formats"></a>134 | [One polyhedral mesh through five formats](gallery/meshes.md#five-formats) | [Meshes](gallery/meshes.md) |
| <a id="su2-zones"></a>138 | [Taylor–Couette sliding case from a two-zone SU2 file](gallery/meshes.md#su2-zones) | [Meshes](gallery/meshes.md) |
| <a id="outlet-backflow"></a>140 | [Turbulence of the backflow at a pressure outlet: an outlet inside the recirculation](gallery/turbulence.md#outlet-backflow) | [Turbulence models](gallery/turbulence.md) |
| <a id="k-sink"></a>141 | [The sink of k near walls of the ζ-f models: grid convergence at the wall](gallery/turbulence.md#k-sink) | [Turbulence models](gallery/turbulence.md) |
| <a id="energy-fluid"></a>154 | [A heated cavity, conduction with k(T) and a Couette flow heated by its friction](gallery/energy.md#energy-fluid) | [Energy equation and heat transfer](gallery/energy.md) |
| <a id="cht-solids"></a>156 | [A heated wall cooled by a channel flow, and a composite wall](gallery/energy.md#cht-solids) | [Energy equation and heat transfer](gallery/energy.md) |
| <a id="cht-turbulent"></a>158 | [A finned heat sink in a turbulent air channel, and a turbulent channel over a conducting wall](gallery/energy.md#cht-turbulent) | [Energy equation and heat transfer](gallery/energy.md) |
| <a id="cht-transient"></a>159 | [A heated wall switched on, and two bodies brought into contact](gallery/energy.md#cht-transient) | [Energy equation and heat transfer](gallery/energy.md) |
| <a id="cht-contact"></a>160 | [A chip on a circuit board with a thermal pad, and an orthotropic patch test](gallery/energy.md#cht-contact) | [Energy equation and heat transfer](gallery/energy.md) |
| <a id="budgets"></a>166 | [Iteration budget, diagnosis and study budgets](gallery/runs.md#budgets) | [Run control, health and adjoint readiness](gallery/runs.md) |
| <a id="limitations"></a>170 | [Verification grades and the limitations of a run: the raw Gmsh pipe](gallery/verification.md#limitations) | [Verification](gallery/verification.md) |
| <a id="tool-registry"></a>171 | [Every function as a tool: the catalogue an agent reads](gallery/runs.md#tool-registry) | [Run control, health and adjoint readiness](gallery/runs.md) |
| <a id="journal"></a>172 | [A session that can be replayed: the journal of tool calls](gallery/runs.md#journal) | [Run control, health and adjoint readiness](gallery/runs.md) |
| <a id="mcp"></a>175 | [An assistant at the controls: the MCP server](gallery/runs.md#mcp) | [Run control, health and adjoint readiness](gallery/runs.md) |
| <a id="control-port"></a>177 | [An assistant steers a running case: the control port](gallery/runs.md#control-port) | [Run control, health and adjoint readiness](gallery/runs.md) |
| <a id="cht-acceleration"></a>180 | [A heated wall brought to its steady and periodic state with thermal acceleration](gallery/energy.md#cht-acceleration) | [Energy equation and heat transfer](gallery/energy.md) |
| <a id="heat-sink-study"></a>181 | [An assistant runs a study: the fin pitch of a heat sink](gallery/runs.md#heat-sink-study) | [Run control, health and adjoint readiness](gallery/runs.md) |
| <a id="mms-turbulence"></a>184 | [Manufactured solutions of k-ε-ζ-f₀, the curvature correction and the QCR: orders and two findings](gallery/turbulence.md#mms-turbulence) | [Turbulence models](gallery/turbulence.md) |
| <a id="buoyancy"></a>186 | [A warm strip in a closed air enclosure, and a stratified box that stays at rest](gallery/energy.md#buoyancy) | [Energy equation and heat transfer](gallery/energy.md) |
| <a id="operator-check"></a>189 | [Which cells can the curvature correction and the QCR trust? The operator check](gallery/meshes.md#operator-check) | [Meshes](gallery/meshes.md) |
| <a id="cht-cycle"></a>190 | [A pulsating channel over a heated wall, coupled by the cycle-averaged loop](gallery/energy.md#cht-cycle) | [Energy equation and heat transfer](gallery/energy.md) |
| <a id="pipe-turbulent"></a>193 | [Turbulent pipe flow at Re_b 5300: friction against the DNS, three meshes with GCI](gallery/turbulence.md#pipe-turbulent) | [Turbulence models](gallery/turbulence.md) |
| <a id="qcr-consistency"></a>194 | [The QCR stress on tetrahedra and perturbed meshes: a consistent divergence](gallery/turbulence.md#qcr-consistency) | [Turbulence models](gallery/turbulence.md) |
| <a id="variable-density"></a>195 | [Air between walls at 960 K and 240 K, a heated duct and two closed boxes](gallery/energy.md#variable-density) | [Energy equation and heat transfer](gallery/energy.md) |
| <a id="tjunction-mixing"></a>196 | [T-junction mixing on polyhedra: three meshes with grid convergence](gallery/flow.md#tjunction-mixing) | [Incompressible flow](gallery/flow.md) |
| <a id="curvature-consistency"></a>197 | [The curvature correction on tetrahedra and perturbed meshes: a consistent DS/Dt](gallery/turbulence.md#curvature-consistency) | [Turbulence models](gallery/turbulence.md) |
| <a id="diffuser"></a>199 | [Asymmetric plane diffuser: separation on one wall, rounded corners and three meshes with GCI](gallery/turbulence.md#diffuser) | [Turbulence models](gallery/turbulence.md) |
| <a id="canonical-ranks"></a>200 | [Three ranks, the bits of one: the canonical local mesh](gallery/meshes.md#canonical-ranks) | [Meshes](gallery/meshes.md) |
| <a id="cht-cycle-skin"></a>202 | [A thin-skin wall under a pulsating channel: the skin swing, fully coupled cycles at the end and the CHT restart](gallery/energy.md#cht-cycle-skin) | [Energy equation and heat transfer](gallery/energy.md) |
| <a id="ubend"></a>204 | [U-bend of Rc/H = 1: separation behind a 180° turn, three meshes with GCI, sharp edges](gallery/turbulence.md#ubend) | [Turbulence models](gallery/turbulence.md) |
| <a id="two-fluids"></a>206 | [Water and oil in counterflow across a steel plate: two fluids in one conjugate case](gallery/energy.md#two-fluids) | [Energy equation and heat transfer](gallery/energy.md) |
| <a id="cht-frozen-flow"></a>208 | [A heated wall warmed up with the flow held, and a pulsating channel warmed up by the cycle loop](gallery/energy.md#cht-frozen-flow) | [Energy equation and heat transfer](gallery/energy.md) |
| <a id="mlp-limiter"></a>211 | [A step that stays a step: the MLP vertex limiter](gallery/scalars.md#mlp-limiter) | [Passive scalars and the wall distance](gallery/scalars.md) |
| <a id="exact-reductions"></a>219 | [Sums that do not care how many ranks there are: exact reductions](gallery/meshes.md#exact-reductions) | [Meshes](gallery/meshes.md) |
| <a id="absolute-pressure"></a>220 | [Absolute pressures: two fluids at 1 bar and 1.5 bar](gallery/boundary-conditions.md#absolute-pressure) | [Boundary conditions](gallery/boundary-conditions.md) |
| <a id="start-values"></a>221 | [Start values per region and per fluid part: every domain of an exchanger starts on its own](gallery/energy.md#start-values) | [Energy equation and heat transfer](gallery/energy.md) |
| <a id="curved-blades"></a>225 | [Three stirrers, one with seven curved blades: circular-arc blades of constant thickness from Gmsh](gallery/flow.md#curved-blades) | [Incompressible flow](gallery/flow.md) |
| <a id="scalar-ranks"></a>230 | [Five ranks, one solution: the scalar solver on ranks](gallery/meshes.md#scalar-ranks) | [Meshes](gallery/meshes.md) |
| <a id="multigrid-ranks"></a>231 | [Five ranks, one hierarchy: multigrid on ranks](gallery/meshes.md#multigrid-ranks) | [Meshes](gallery/meshes.md) |
| <a id="flow-ranks"></a>234 | [Five ranks, one flow: SIMPLEC on ranks](gallery/meshes.md#flow-ranks) | [Meshes](gallery/meshes.md) |
| <a id="ledger"></a>238 | [The conservation ledger as data: a heated cavity's mass, kinetic and thermal accounts](gallery/energy.md#ledger) | [Energy equation and heat transfer](gallery/energy.md) |
| <a id="spray-duct"></a>242 | [Water droplets sprayed into an air duct](gallery/sprays.md#spray-duct) | [Lagrangian sprays](gallery/sprays.md) |
| <a id="flow-models-ranks"></a>244 | [Five ranks, one turbulent heated flow: models on ranks](gallery/meshes.md#flow-models-ranks) | [Meshes](gallery/meshes.md) |
| <a id="furnace"></a>259 | [A box furnace with cooled tubes, and the cases with exact solutions](gallery/radiation.md#furnace) | [Thermal radiation](gallery/radiation.md) |
| <a id="restart-ranks"></a>262 | [Written on five ranks, read on any: restart files on ranks](gallery/meshes.md#restart-ranks) | [Meshes](gallery/meshes.md) |
| <a id="spray-datums"></a>265 | [Where the liquid goes: datums and droplet statistics](gallery/sprays.md#spray-datums) | [Lagrangian sprays](gallery/sprays.md) |
| <a id="penalised-cylinder"></a>270 | [A cylinder put into a box mesh by Brinkman penalisation](gallery/immersed.md#penalised-cylinder) | [Immersed bodies](gallery/immersed.md) |
| <a id="adjoint-a0"></a>271 | [Adjoint stage A0: one steady residual for SIMPLEC and the coupled solver](gallery/verification.md#adjoint-a0) | [Verification](gallery/verification.md) |
| <a id="natural-convection"></a>272 | [Natural convection benchmarks: a heated square cavity at Ra 10³ to 10⁶, a vertical plate and mixed convection with flow reversal](gallery/energy.md#natural-convection) | [Energy equation and heat transfer](gallery/energy.md) |
| <a id="porous-flow"></a>273 | [A ceramic honeycomb in an air channel, and a porous block against the exact pressure drop](gallery/porous.md#porous-flow) | [Porous regions](gallery/porous.md) |
| <a id="porous-two-temperatures"></a>275 | [Charging a honeycomb heat store, and Schumann's bed against the exact solution](gallery/porous.md#porous-two-temperatures) | [Porous regions](gallery/porous.md) |
| <a id="species-transport"></a>278 | [A methane jet in a coflow of air, and binary diffusion against the exact solution](gallery/species.md#species-transport) | [Species and reacting flows](gallery/species.md) |
| <a id="run-ranks"></a>279 | [Solved on five ranks, written by one: a run on ranks](gallery/meshes.md#run-ranks) | [Meshes](gallery/meshes.md) |
| <a id="spray-two-way"></a>280 | [A water spray that drives the air](gallery/sprays.md#spray-two-way) | [Lagrangian sprays](gallery/sprays.md) |
| <a id="furnace-paths"></a>281 | [Paths of energy bundles in the box furnace, and the net flux into every wall](gallery/radiation.md#furnace-paths) | [Thermal radiation](gallery/radiation.md) |
| <a id="spray-cooling"></a>283 | [Water droplets cool a hot air stream](gallery/sprays.md#spray-cooling) | [Lagrangian sprays](gallery/sprays.md) |
| <a id="spray-dispersion"></a>285 | [Droplets spread by the turbulence of a duct](gallery/sprays.md#spray-dispersion) | [Lagrangian sprays](gallery/sprays.md) |
| <a id="gas-chemistry"></a>289 | [Methane ignites in a constant-pressure reactor and ends at its chemical equilibrium](gallery/species.md#gas-chemistry) | [Species and reacting flows](gallery/species.md) |
| <a id="run-ranks-physics"></a>290 | [Five ranks, one turbulent heated transient: the physics of a run on ranks](gallery/meshes.md#run-ranks-physics) | [Meshes](gallery/meshes.md) |
| <a id="topology-bend"></a>291 | [Topology search of a pipe bend, and how accurate the penalised wall is](gallery/immersed.md#topology-bend) | [Immersed bodies](gallery/immersed.md) |
| <a id="spray-steady"></a>292 | [The spray-cooled duct as a steady state](gallery/sprays.md#spray-steady) | [Lagrangian sprays](gallery/sprays.md) |
| <a id="reed-beam"></a>295 | [A reed valve as a beam: modes, and the clamp spring against the frequency equation](gallery/structure.md#reed-beam) | [Structures](gallery/structure.md) |
| <a id="deforming-mesh"></a>296 | [A piston compresses a closed gas volume on a deforming mesh](gallery/moving-solids.md#deforming-mesh) | [Moving meshes and solids](gallery/moving-solids.md) |
| <a id="reed-plate"></a>297 | [A reed petal as a shell: bending and twist on clamp springs](gallery/structure.md#reed-plate) | [Structures](gallery/structure.md) |
| <a id="reed-dynamics"></a>299 | [A reed struck by a pressure pulse, and a hot reed](gallery/structure.md#reed-dynamics) | [Structures](gallery/structure.md) |
| <a id="run-ranks-mrf"></a>301 | [A stirred vessel on five ranks: rotating zones of MRF on ranks](gallery/meshes.md#run-ranks-mrf) | [Meshes](gallery/meshes.md) |

## Adding a feature

Every feature ships with a gallery example, case files and a recipe (D18,
as in PyAeroFlow); see [CONTRIBUTING.md](../CONTRIBUTING.md). Put the example on the topic page that fits, or start a new page when a new physics or topic begins (energy, turbulence, reacting flows, sprays, ...). Number it with its pull request (a second example of the same pull request gets a letter: 54b), give it an `<a id>` anchor and a row in that page's table of contents. The two tables above are regenerated after the merge (`python -m tools.docs.indexes`); pull requests leave them alone. A new topic page also needs a row in the first table.
