# Recipes

Short command sequences for everyday tasks. Each recipe names the case-file
keys it uses; the keys of the benchmark cases are listed in the
[Case catalogue](cases.md#case-file-keys). Pictures of complete cases are in
the [Example gallery](gallery.md).

The recipes are split by topic; each page has its own table of contents. The numbers are global and stay the same, so "recipe 12" always means the same entry.

| Page | What it covers | Recipes |
|---|---|---|
| [Kernel benchmark G1](recipes/benchmarks.md) | build and run gate G1 on a CPU, on Windows and on an NVIDIA GPU, and rerun it on unstructured meshes | 1, 2, 3, 10 |
| [Development: tests, CI and releases](recipes/development.md) | MPI tests, checking a change the way CI does, releases, derivative checks of the kernels, manufactured solutions of the turbulence models | 4, 5, 6, 25, 184 |
| [Meshes](recipes/meshes.md) | generate, read, convert and partition meshes, update the mesher's golden values, run flows on poor Gmsh meshes, the operator exactness check, the canonical local mesh and exact reductions | 7, 8, 9, 11, 36, 127, 129, 134, 138, 189, 200, 219, 230, 231, 234, 244, 262, 279, 290, 301, 304, 305 |
| [Scalar transport and the pressure solver](recipes/scalar-pressure.md) | order of accuracy of the scalar operators, convection schemes, pressure-solver studies, hypre | 12, 13, 15, 17 |
| [GPU runs](recipes/gpu.md) | scalar operators, pressure solver and flow cases on an NVIDIA GPU (also a free Colab T4) | 14, 16, 24 |
| [Incompressible flow](recipes/flow.md) | Couette, Poiseuille, cavity and cylinder runs, boundary conditions, transient runs, time-step control, the coupled algorithm, heat transfer in turbulent flow, rotating machinery, the T-junction mixing study | 18, 19, 20, 21, 23, 28, 78, 80, 81, 83, 85, 86, 87, 94, 95, 96, 98, 99, 100, 101, 102, 225, 103, 196 |
| [Passive scalars and the wall distance](recipes/scalars.md) | passive scalars with their boundary conditions, bounds, exact solutions, the Nusselt number of a heated channel, the wall distance | 63, 211 |
| [Energy equation and heat transfer](recipes/energy.md) | the temperature of a flow: properties of T, thermal wall conditions, viscous heating, the heat balance, solid regions and conjugate heat transfer, laminar, turbulent and transient, thermal acceleration, the cycle loop and the frozen-flow march, several fluids, start values per region, buoyancy and variable density | 154, 156, 158, 159, 160, 180, 190, 202, 208, 206, 186, 195, 272, 221, 238 |
| [Thermal radiation](recipes/radiation.md) | the radiation of a prescribed temperature field by Monte Carlo bundles: walls, noise, threads, ranks and the GPU, the checks against exact solutions | 259, 281 |
| [Turbulence models](recipes/turbulence.md) | RANS turbulence models: the model, the wall treatment by the mesh, inlet turbulence, roughness, the near-wall checks, curvature correction and QCR on any mesh, the internal flows of D45 (pipe, diffuser, U-bend) | 68, 71, 73, 74, 88, 89, 90, 113, 115, 141, 140, 194, 193, 199, 197, 204 |
| [Boundary conditions](recipes/boundary-conditions.md) | slip walls, wedges for axisymmetric flow, inlet profiles from tables, backflow at pressure outlets, convective outlets, periodic modules, absolute or gauge pressure values (D239) | 29, 67, 220 |
| [Lagrangian sprays](recipes/sprays.md) | spray droplets into a flow with an injector file, wall behaviour, the mass book | 242, 265, 280, 283, 285, 292, 306 |
| [Porous regions](recipes/porous.md) | porous blocks and honeycombs in a flow: resistances, axes, the honeycomb shorthand | 273, 275 |
| [Two fluids (volume of fluid)](recipes/multiphase.md) | water and air in one case: phases, start regions, inflow, `vof.csv` | 317 |
| [Structures](recipes/structure.md) | a structural case of beams, plates, solids and clamp springs: static and modal analyses, the clamp spring fitted to a measured frequency, runs in time with damping and temperatures, contact with a seat and a stop, temperatures from a flow case, the results | 295, 297, 299, 300, 302, 308 |
| [Immersed bodies](recipes/immersed.md) | put a body in the flow without meshing it: `[immersed]`, the permeability, forces | 270, 291, 310 |
| [Moving meshes and solids](recipes/moving-solids.md) | move the mesh with a prescribed displacement: `[mesh_motion]`, moving walls, a piston on a closed volume | 296 |
| [Species and reacting flows](recipes/species.md) | species of an ideal-gas mixture: names, compositions in mole or mass fractions, diffusion by Schmidt or Lewis numbers, mixture properties, the element ledger | 278, 289, 298 |
| [Expressions and materials](recipes/properties.md) | values as expressions of x, y, z and t, constants and tables, wall mass flux, the material database, regions | 31, 32, 277 |
| [Python](recipes/python.md) | the `foehncraft` package: cases, expressions, runs with callbacks, results, materials, studies | 33 |
| [Verification](recipes/verification.md) | the discretisation uncertainty of a sweep (GCI), the generated verification manual | 34, 35, 165, 170, 271 |
| [Post-processing](recipes/post-processing.md) | derived fields with the output and from restart files (`foehn post derive`) | 30, 108, 110 |
| [Engineering gallery](recipes/applications.md) | run an application of the Engineering Gallery and compare it with published data | 161, 203 |
| [Run control and health](recipes/runs.md) | restart, stop and steer runs, drive runs from a script or an agent, balances and divergence reports, the log, scheduled changes and case upgrades, iteration budgets, the tool registry, journals, the MCP server and the control port | 22, 26, 27, 37, 38, 39, 166, 171, 172, 175, 177, 181, 316 |

## All recipes by number

| # | Recipe | Page |
|---|---|---|
| <a id="g1-cpu"></a>1 | [Build, test and run the G1 benchmark on your CPU](recipes/benchmarks.md#g1-cpu) | [Kernel benchmark G1](recipes/benchmarks.md) |
| <a id="g1-windows"></a>2 | [Run G1 on Windows: CPU verdict and the CUDA reference](recipes/benchmarks.md#g1-windows) | [Kernel benchmark G1](recipes/benchmarks.md) |
| <a id="g1-gpu"></a>3 | [Run gate G1 on an NVIDIA GPU](recipes/benchmarks.md#g1-gpu) | [Kernel benchmark G1](recipes/benchmarks.md) |
| <a id="mpi"></a>4 | [Run the MPI test with several ranks on one machine](recipes/development.md#mpi) | [Development: tests, CI and releases](recipes/development.md) |
| <a id="ci-local"></a>5 | [Check a change the way CI does](recipes/development.md#ci-local) | [Development: tests, CI and releases](recipes/development.md) |
| <a id="release"></a>6 | [Make a release (tag and GitHub release)](recipes/development.md#release) | [Development: tests, CI and releases](recipes/development.md) |
| <a id="mesh-generate"></a>7 | [Generate a test mesh and check its quality](recipes/meshes.md#mesh-generate) | [Meshes](recipes/meshes.md) |
| <a id="mesh-gmsh"></a>8 | [Read a Gmsh mesh and convert it into polyhedra](recipes/meshes.md#mesh-gmsh) | [Meshes](recipes/meshes.md) |
| <a id="mesh-partition"></a>9 | [Partition a mesh and test the halo exchange](recipes/meshes.md#mesh-partition) | [Meshes](recipes/meshes.md) |
| <a id="g1-stage1"></a>10 | [Rerun G1 on polyhedral and tetrahedral meshes](recipes/benchmarks.md#g1-stage1) | [Kernel benchmark G1](recipes/benchmarks.md) |
| <a id="golden"></a>11 | [Update the golden values of the mesher](recipes/meshes.md#golden) | [Meshes](recipes/meshes.md) |
| <a id="scalar-mms"></a>12 | [Check the order of accuracy of the scalar operators](recipes/scalar-pressure.md#scalar-mms) | [Scalar transport and the pressure solver](recipes/scalar-pressure.md) |
| <a id="scalar-schemes"></a>13 | [Compare convection schemes and limiters](recipes/scalar-pressure.md#scalar-schemes) | [Scalar transport and the pressure solver](recipes/scalar-pressure.md) |
| <a id="scalar-gpu"></a>14 | [Run the scalar operators on an NVIDIA GPU](recipes/gpu.md#scalar-gpu) | [GPU runs](recipes/gpu.md) |
| <a id="poisson-study"></a>15 | [Study the pressure solver on a Poisson problem](recipes/scalar-pressure.md#poisson-study) | [Scalar transport and the pressure solver](recipes/scalar-pressure.md) |
| <a id="poisson-gpu"></a>16 | [Run the pressure solver on an NVIDIA GPU against hypre](recipes/gpu.md#poisson-gpu) | [GPU runs](recipes/gpu.md) |
| <a id="hypre"></a>17 | [Build hypre and compare with it](recipes/scalar-pressure.md#hypre) | [Scalar transport and the pressure solver](recipes/scalar-pressure.md) |
| <a id="flow-verification"></a>18 | [Verify the flow solver on Couette and Poiseuille flow](recipes/flow.md#flow-verification) | [Incompressible flow](recipes/flow.md) |
| <a id="flow-cavity"></a>19 | [Run the lid-driven cavity and compare with Ghia et al.](recipes/flow.md#flow-cavity) | [Incompressible flow](recipes/flow.md) |
| <a id="flow-boundaries"></a>20 | [Set up the boundary conditions of a flow case](recipes/flow.md#flow-boundaries) | [Incompressible flow](recipes/flow.md) |
| <a id="transient-flow"></a>21 | [Run a transient flow case and check its time order](recipes/flow.md#transient-flow) | [Incompressible flow](recipes/flow.md) |
| <a id="restart"></a>22 | [Restart a run, stop it and steer it while it runs](recipes/runs.md#restart) | [Run control and health](recipes/runs.md) |
| <a id="cylinder"></a>23 | [Run the cylinder benchmark of Schäfer & Turek](recipes/flow.md#cylinder) | [Incompressible flow](recipes/flow.md) |
| <a id="flow-gpu"></a>24 | [Run a flow case on an NVIDIA GPU and compare it with the CPU](recipes/gpu.md#flow-gpu) | [GPU runs](recipes/gpu.md) |
| <a id="adjoint-ready"></a>25 | [Check the derivatives of the kernels (adjoint readiness)](recipes/development.md#adjoint-ready) | [Development: tests, CI and releases](recipes/development.md) |
| <a id="run-contract"></a>26 | [Drive runs from a script or an agent](recipes/runs.md#run-contract) | [Run control and health](recipes/runs.md) |
| <a id="health"></a>27 | [Check the balances and read a divergence report](recipes/runs.md#health) | [Run control and health](recipes/runs.md) |
| <a id="time-control"></a>28 | [Let the error choose the time step, write output at exact times, start from potential flow](recipes/flow.md#time-control) | [Incompressible flow](recipes/flow.md) |
| <a id="boundary-conditions"></a>29 | [Use slip walls, wedges, table profiles, backflow treatment and convective outlets](recipes/boundary-conditions.md#boundary-conditions) | [Boundary conditions](recipes/boundary-conditions.md) |
| <a id="derived-fields"></a>30 | [Write vorticity, Q, λ₂, wall shear stress and other derived fields](recipes/post-processing.md#derived-fields) | [Post-processing](recipes/post-processing.md) |
| <a id="expressions"></a>31 | [Write boundary values, body forces and the viscosity as expressions](recipes/properties.md#expressions) | [Expressions and materials](recipes/properties.md) |
| <a id="materials"></a>32 | [Take the fluid properties from the material database](recipes/properties.md#materials) | [Expressions and materials](recipes/properties.md) |
| <a id="python-api"></a>33 | [Build, run and study cases from Python](recipes/python.md#python-api) | [Python](recipes/python.md) |
| <a id="gci"></a>34 | [Estimate the discretisation uncertainty of a mesh or time-step sweep](recipes/verification.md#gci) | [Verification](recipes/verification.md) |
| <a id="manual"></a>35 | [Regenerate the verification manual](recipes/verification.md#manual) | [Verification](recipes/verification.md) |
| <a id="skewed-tets"></a>36 | [Run a flow on raw Gmsh tetrahedra (skewed faces)](recipes/meshes.md#skewed-tets) | [Meshes](recipes/meshes.md) |
| <a id="log"></a>37 | [Read the log and upgrade a case file](recipes/runs.md#log) | [Run control and health](recipes/runs.md) |
| <a id="log-columns"></a>38 | [Choose the columns of the log table](recipes/runs.md#log-columns) | [Run control and health](recipes/runs.md) |
| <a id="schedule"></a>39 | [Prepare changes of settings before the run](recipes/runs.md#schedule) | [Run control and health](recipes/runs.md) |
| <a id="passive-scalars"></a>63 | [Carry a passive scalar and compute a Nusselt number](recipes/scalars.md#passive-scalars) | [Passive scalars and the wall distance](recipes/scalars.md) |
| <a id="periodic"></a>67 | [Compute a fully developed flow on one periodic module](recipes/boundary-conditions.md#periodic) | [Boundary conditions](recipes/boundary-conditions.md) |
| <a id="sst"></a>68 | [Run a case with the k-ω SST model](recipes/turbulence.md#sst) | [Turbulence models](recipes/turbulence.md) |
| <a id="k-epsilon"></a>71 | [Run a case with the k-ε family and wall functions](recipes/turbulence.md#k-epsilon) | [Turbulence models](recipes/turbulence.md) |
| <a id="zeta-f"></a>73 | [Run a case with k-ε-ζ-f and resolved walls](recipes/turbulence.md#zeta-f) | [Turbulence models](recipes/turbulence.md) |
| <a id="zeta-f0"></a>74 | [Run a case with k-ε-ζ-f₀ at any first-cell y⁺](recipes/turbulence.md#zeta-f0) | [Turbulence models](recipes/turbulence.md) |
| <a id="coupled"></a>78 | [Converge a steady case in fewer iterations with the coupled algorithm](recipes/flow.md#coupled) | [Incompressible flow](recipes/flow.md) |
| <a id="momentum-multigrid"></a>80 | [Choose and tune the momentum sweep of the coupled algorithm](recipes/flow.md#momentum-multigrid) | [Incompressible flow](recipes/flow.md) |
| <a id="turbulence-multigrid"></a>81 | [Use the multigrid of the turbulence equations and its larger steps](recipes/flow.md#turbulence-multigrid) | [Incompressible flow](recipes/flow.md) |
| <a id="pressure-hierarchy"></a>83 | [See and tune the pressure hierarchy of the coupled preconditioner](recipes/flow.md#pressure-hierarchy) | [Incompressible flow](recipes/flow.md) |
| <a id="thermal-wall-function"></a>85 | [Add heat transfer to a turbulent flow: wall function and turbulent Prandtl number](recipes/flow.md#thermal-wall-function) | [Incompressible flow](recipes/flow.md) |
| <a id="heated-channel-dns"></a>86 | [Compare the heated channel with the DNS Nusselt numbers](recipes/flow.md#heated-channel-dns) | [Incompressible flow](recipes/flow.md) |
| <a id="channel-dns-reynolds"></a>87 | [Run the heated channel at another Re_τ or Pr](recipes/flow.md#channel-dns-reynolds) | [Incompressible flow](recipes/flow.md) |
| <a id="elliptic-blending"></a>88 | [Run a case with the elliptic-blending k-ε](recipes/turbulence.md#elliptic-blending) | [Turbulence models](recipes/turbulence.md) |
| <a id="step-validation"></a>89 | [Validate the turbulence models on the backward-facing step](recipes/turbulence.md#step-validation) | [Turbulence models](recipes/turbulence.md) |
| <a id="impinging-jet"></a>90 | [Feed a jet with a developed pipe profile and compute the Nusselt number](recipes/turbulence.md#impinging-jet) | [Turbulence models](recipes/turbulence.md) |
| <a id="mrf"></a>94 | [Run a rotating zone in a fixed frame (pump, fan, mixer, rotor)](recipes/flow.md#mrf) | [Incompressible flow](recipes/flow.md) |
| <a id="sliding-mesh"></a>95 | [Turn a rotor through the mesh in time (sliding mesh, conformal stepping)](recipes/flow.md#sliding-mesh) | [Incompressible flow](recipes/flow.md) |
| <a id="sliding-mesh-imprint"></a>96 | [Turn a rotor against a stator with a different pitch (sliding mesh by imprinting)](recipes/flow.md#sliding-mesh-imprint) | [Incompressible flow](recipes/flow.md) |
| <a id="stirred-vessel"></a>98 | [Mesh and run a rotor with blades from Gmsh, and draw its relative flow](recipes/flow.md#stirred-vessel) | [Incompressible flow](recipes/flow.md) |
| <a id="rotor-fields"></a>99 | [Write the relative flow of a rotor: relative and frame velocity, flow angles, rothalpy](recipes/flow.md#rotor-fields) | [Incompressible flow](recipes/flow.md) |
| <a id="turbulent-rotor"></a>100 | [Run a turbulent rotor: the frame or the sliding mesh with a turbulence model](recipes/flow.md#turbulent-rotor) | [Incompressible flow](recipes/flow.md) |
| <a id="torque-power"></a>101 | [Read the torque and power of a rotor, and the head and efficiency of a machine](recipes/flow.md#torque-power) | [Incompressible flow](recipes/flow.md) |
| <a id="several-rotors"></a>102 | [Run several rotors at their own speeds, and a transient frozen rotor](recipes/flow.md#several-rotors) | [Incompressible flow](recipes/flow.md) |
| <a id="revolved-interfaces"></a>103 | [Give a rotor an interface of any rotationally symmetric shape: discs, cones, pots](recipes/flow.md#revolved-interfaces) | [Incompressible flow](recipes/flow.md) |
| <a id="turbulence-derived-fields"></a>108 | [Write the Reynolds stress, the scales and the wall values of a turbulence model](recipes/post-processing.md#turbulence-derived-fields) | [Post-processing](recipes/post-processing.md) |
| <a id="transient-statistics"></a>110 | [Average a transient run: mean, RMS, covariance, confidence interval](recipes/post-processing.md#transient-statistics) | [Post-processing](recipes/post-processing.md) |
| <a id="curvature-correction"></a>113 | [Make SST sensitive to curvature and system rotation](recipes/turbulence.md#curvature-correction) | [Turbulence models](recipes/turbulence.md) |
| <a id="qcr"></a>115 | [Get the secondary flow of non-circular ducts with SST](recipes/turbulence.md#qcr) | [Turbulence models](recipes/turbulence.md) |
| <a id="mesh-su2"></a>127 | [Read an SU2 mesh and map its markers to boundary conditions](recipes/meshes.md#mesh-su2) | [Meshes](recipes/meshes.md) |
| <a id="mesh-cgns-openfoam"></a>129 | [Read a CGNS or OpenFOAM mesh](recipes/meshes.md#mesh-cgns-openfoam) | [Meshes](recipes/meshes.md) |
| <a id="mesh-msh-vtk"></a>134 | [Read a Fluent `.msh`, VTK or VTKHDF mesh, write a VTU](recipes/meshes.md#mesh-msh-vtk) | [Meshes](recipes/meshes.md) |
| <a id="mesh-multizone"></a>138 | [Run a rotor/stator sliding case from a multi-zone SU2 or CGNS file](recipes/meshes.md#mesh-multizone) | [Meshes](recipes/meshes.md) |
| <a id="outlet-backflow"></a>140 | [Say what turbulence enters through an outlet](recipes/turbulence.md#outlet-backflow) | [Turbulence models](recipes/turbulence.md) |
| <a id="k-sink"></a>141 | [Choose the sink of k of the ζ-f models](recipes/turbulence.md#k-sink) | [Turbulence models](recipes/turbulence.md) |
| <a id="energy-fluid"></a>154 | [Solve the temperature of a flow with temperature-dependent properties](recipes/energy.md#energy-fluid) | [Energy equation and heat transfer](recipes/energy.md) |
| <a id="cht-solids"></a>156 | [Couple the temperature of a fluid and a solid (conjugate heat transfer)](recipes/energy.md#cht-solids) | [Energy equation and heat transfer](recipes/energy.md) |
| <a id="cht-turbulent"></a>158 | [Couple a turbulent flow with a solid (turbulent conjugate heat transfer)](recipes/energy.md#cht-turbulent) | [Energy equation and heat transfer](recipes/energy.md) |
| <a id="cht-transient"></a>159 | [Follow how a solid and a fluid heat up together (transient conjugate heat transfer)](recipes/energy.md#cht-transient) | [Energy equation and heat transfer](recipes/energy.md) |
| <a id="cht-contact"></a>160 | [Put a contact resistance between two parts, orient an orthotropic solid, heat a component](recipes/energy.md#cht-contact) | [Energy equation and heat transfer](recipes/energy.md) |
| <a id="fyfe-stirred-tank-recipe"></a>161 | [Run the 2 L stirred tank of Fyfe et al.: power number and mixing time](recipes/applications.md#fyfe-stirred-tank-recipe) | [Engineering gallery](recipes/applications.md) |
| <a id="status"></a>165 | [Grade a feature on the verification status page](recipes/verification.md#status) | [Verification](recipes/verification.md) |
| <a id="budgets"></a>166 | [Set an iteration budget and a budget per study](recipes/runs.md#budgets) | [Run control and health](recipes/runs.md) |
| <a id="limitations"></a>170 | [Read and act on the limitations of a run](recipes/verification.md#limitations) | [Verification](recipes/verification.md) |
| <a id="tools"></a>171 | [Call any tool from a script or an agent](recipes/runs.md#tools) | [Run control and health](recipes/runs.md) |
| <a id="journal"></a>172 | [Record a session and replay it](recipes/runs.md#journal) | [Run control and health](recipes/runs.md) |
| <a id="mcp"></a>175 | [Connect an assistant through the MCP server](recipes/runs.md#mcp) | [Run control and health](recipes/runs.md) |
| <a id="control-port"></a>177 | [Watch and steer a running case](recipes/runs.md#control-port) | [Run control and health](recipes/runs.md) |
| <a id="cht-acceleration"></a>180 | [Bring the solids of a transient conjugate run to their state sooner (thermal acceleration)](recipes/energy.md#cht-acceleration) | [Energy equation and heat transfer](recipes/energy.md) |
| <a id="heat-sink-study"></a>181 | [Let an assistant run a parameter study](recipes/runs.md#heat-sink-study) | [Run control and health](recipes/runs.md) |
| <a id="mms-turbulence"></a>184 | [Check a turbulence model with a manufactured solution](recipes/development.md#mms-turbulence) | [Development: tests, CI and releases](recipes/development.md) |
| <a id="buoyancy"></a>186 | [Switch on gravity and buoyancy; which pressure do my boundary values mean](recipes/energy.md#buoyancy) | [Energy equation and heat transfer](recipes/energy.md) |
| <a id="operator-check"></a>189 | [Check which cells the velocity gradient and Hessian are exact on](recipes/meshes.md#operator-check) | [Meshes](recipes/meshes.md) |
| <a id="cht-cycle"></a>190 | [Couple a periodic flow with slow solids by cycle averages (exhaust manifolds, pulsating flows)](recipes/energy.md#cht-cycle) | [Energy equation and heat transfer](recipes/energy.md) |
| <a id="pipe-turbulent"></a>193 | [Run a turbulent pipe and check its friction](recipes/turbulence.md#pipe-turbulent) | [Turbulence models](recipes/turbulence.md) |
| <a id="qcr-consistency"></a>194 | [Run the QCR on tetrahedra and perturbed meshes](recipes/turbulence.md#qcr-consistency) | [Turbulence models](recipes/turbulence.md) |
| <a id="variable-density"></a>195 | [Let the density of a gas follow its temperature (low-Mach variable density)](recipes/energy.md#variable-density) | [Energy equation and heat transfer](recipes/energy.md) |
| <a id="tjunction-mixing"></a>196 | [Mix two streams in a T-junction and check the mesh](recipes/flow.md#tjunction-mixing) | [Incompressible flow](recipes/flow.md) |
| <a id="curvature-consistency"></a>197 | [Run the curvature correction on tetrahedra and perturbed meshes](recipes/turbulence.md#curvature-consistency) | [Turbulence models](recipes/turbulence.md) |
| <a id="diffuser"></a>199 | [Run an asymmetric diffuser, its grid study and a transient run from the steady state](recipes/turbulence.md#diffuser) | [Turbulence models](recipes/turbulence.md) |
| <a id="canonical-ranks"></a>200 | [Check that every rank has the bits of one rank (canonical local mesh)](recipes/meshes.md#canonical-ranks) | [Meshes](recipes/meshes.md) |
| <a id="cht-cycle-skin"></a>202 | [Get the swing of the wall temperature from the cycle loop, and restart a conjugate run](recipes/energy.md#cht-cycle-skin) | [Energy equation and heat transfer](recipes/energy.md) |
| <a id="shell-tube-hx-recipe"></a>203 | [Mesh and run a multi-region heat exchanger from a STEP assembly](recipes/applications.md#shell-tube-hx-recipe) | [Engineering gallery](recipes/applications.md) |
| <a id="ubend"></a>204 | [Run a U-bend, its grid study and the sharp-edged variant](recipes/turbulence.md#ubend) | [Turbulence models](recipes/turbulence.md) |
| <a id="two-fluids"></a>206 | [Put several fluids in one conjugate case (water and air in a heat exchanger)](recipes/energy.md#two-fluids) | [Energy equation and heat transfer](recipes/energy.md) |
| <a id="cht-frozen-flow"></a>208 | [Warm up a conjugate run with the flow held, or the warm-up history of a periodic flow](recipes/energy.md#cht-frozen-flow) | [Energy equation and heat transfer](recipes/energy.md) |
| <a id="mlp-limiter"></a>211 | [Keep a scalar or the temperature within its bounds on tetrahedra and polyhedra](recipes/scalars.md#mlp-limiter) | [Passive scalars and the wall distance](recipes/scalars.md) |
| <a id="exact-reductions"></a>219 | [Sum over ranks exactly (and check it)](recipes/meshes.md#exact-reductions) | [Meshes](recipes/meshes.md) |
| <a id="absolute-pressure"></a>220 | [Give pressures as absolute values, or keep gauge values](recipes/boundary-conditions.md#absolute-pressure) | [Boundary conditions](recipes/boundary-conditions.md) |
| <a id="start-values"></a>221 | [Start every region and fluid part on its own](recipes/energy.md#start-values) | [Energy equation and heat transfer](recipes/energy.md) |
| <a id="curved-blades"></a>225 | [Mesh a rotor with curved blades of constant thickness in Gmsh](recipes/flow.md#curved-blades) | [Incompressible flow](recipes/flow.md) |
| <a id="scalar-ranks"></a>230 | [Solve a scalar on ranks (and check it against one rank)](recipes/meshes.md#scalar-ranks) | [Meshes](recipes/meshes.md) |
| <a id="multigrid-ranks"></a>231 | [Use the multigrid preconditioner on ranks](recipes/meshes.md#multigrid-ranks) | [Meshes](recipes/meshes.md) |
| <a id="flow-ranks"></a>234 | [Solve a flow on ranks (and check it against one rank)](recipes/meshes.md#flow-ranks) | [Meshes](recipes/meshes.md) |
| <a id="ledger"></a>238 | [Read the conservation ledger of a run as data, and book a module's rows on it](recipes/energy.md#ledger) | [Energy equation and heat transfer](recipes/energy.md) |
| <a id="spray-duct"></a>242 | [Spray droplets into a flow](recipes/sprays.md#spray-duct) | [Lagrangian sprays](recipes/sprays.md) |
| <a id="flow-models-ranks"></a>244 | [Check a turbulent, heated or transient flow on ranks](recipes/meshes.md#flow-models-ranks) | [Meshes](recipes/meshes.md) |
| <a id="radiation-call"></a>259 | [Compute the radiation of a prescribed temperature field, and check it](recipes/radiation.md#radiation-call) | [Thermal radiation](recipes/radiation.md) |
| <a id="restart-ranks"></a>262 | [Write a restart on ranks and continue on another rank count](recipes/meshes.md#restart-ranks) | [Meshes](recipes/meshes.md) |
| <a id="spray-datums"></a>265 | [Measure the liquid and the droplet sizes at datums](recipes/sprays.md#spray-datums) | [Lagrangian sprays](recipes/sprays.md) |
| <a id="penalised-body"></a>270 | [Put a body in the flow without meshing it](recipes/immersed.md#penalised-body) | [Immersed bodies](recipes/immersed.md) |
| <a id="adjoint-a0"></a>271 | [Evaluate the steady residual and check the reverse sweep of the adjoint](recipes/verification.md#adjoint-a0) | [Verification](recipes/verification.md) |
| <a id="natural-convection"></a>272 | [Set up a natural convection case and check it: a closed cavity, an open plate, a channel with a net flow](recipes/energy.md#natural-convection) | [Energy equation and heat transfer](recipes/energy.md) |
| <a id="porous-flow"></a>273 | [Put a porous block or a honeycomb into a flow](recipes/porous.md#porous-flow) | [Porous regions](recipes/porous.md) |
| <a id="porous-two-temperatures"></a>275 | [Give the matrix a temperature of its own: a packed bed or a heat store](recipes/porous.md#porous-two-temperatures) | [Porous regions](recipes/porous.md) |
| <a id="steam-if97"></a>277 | [Use water and steam of IAPWS-IF97: states, saturation lines and the `if97` material form](recipes/properties.md#steam-if97) | [Expressions and materials](recipes/properties.md) |
| <a id="species-transport"></a>278 | [Carry the species of a gas mixture with the flow](recipes/species.md#species-transport) | [Species and reacting flows](recipes/species.md) |
| <a id="run-ranks"></a>279 | [Run a case on ranks](recipes/meshes.md#run-ranks) | [Meshes](recipes/meshes.md) |
| <a id="spray-two-way"></a>280 | [Let the spray push the gas (two-way coupling)](recipes/sprays.md#spray-two-way) | [Lagrangian sprays](recipes/sprays.md) |
| <a id="radiation-paths"></a>281 | [Follow chosen bundles of a radiation call](recipes/radiation.md#radiation-paths) | [Thermal radiation](recipes/radiation.md) |
| <a id="spray-cooling"></a>283 | [Cool a gas with droplets (two-way heat)](recipes/sprays.md#spray-cooling) | [Lagrangian sprays](recipes/sprays.md) |
| <a id="spray-dispersion"></a>285 | [Let the turbulence spread the droplets](recipes/sprays.md#spray-dispersion) | [Lagrangian sprays](recipes/sprays.md) |
| <a id="gas-chemistry"></a>289 | [Run a mechanism in a constant-pressure reactor and find the equilibrium](recipes/species.md#gas-chemistry) | [Species and reacting flows](recipes/species.md) |
| <a id="run-ranks-physics"></a>290 | [Run a turbulent, heated or transient case on ranks](recipes/meshes.md#run-ranks-physics) | [Meshes](recipes/meshes.md) |
| <a id="topology-search"></a>291 | [Search the shape of a duct with a density field](recipes/immersed.md#topology-search) | [Immersed bodies](recipes/immersed.md) |
| <a id="spray-steady"></a>292 | [Solve a spray in a steady flow (steady spray mode)](recipes/sprays.md#spray-steady) | [Lagrangian sprays](recipes/sprays.md) |
| <a id="structure-case"></a>295 | [Run a structural case of beams and clamp springs, and fit the clamp to a measured frequency](recipes/structure.md#structure-case) | [Structures](recipes/structure.md) |
| <a id="deforming-mesh"></a>296 | [Move the mesh with a prescribed displacement (a piston)](recipes/moving-solids.md#deforming-mesh) | [Moving meshes and solids](recipes/moving-solids.md) |
| <a id="structure-plate"></a>297 | [Model a flat reed as a shell clamped along an edge](recipes/structure.md#structure-plate) | [Structures](recipes/structure.md) |
| <a id="reacting-flow-cells"></a>298 | [Let the species react in the cells of a laminar flow](recipes/species.md#reacting-flow-cells) | [Species and reacting flows](recipes/species.md) |
| <a id="structure-transient"></a>299 | [Run a structure in time, with damping, load histories and temperatures](recipes/structure.md#structure-transient) | [Structures](recipes/structure.md) |
| <a id="structure-solid"></a>300 | [Model a thick part as a solid, from Gmsh or as a block](recipes/structure.md#structure-solid) | [Structures](recipes/structure.md) |
| <a id="run-ranks-mrf"></a>301 | [Run a case with rotating zones (MRF) or a closed domain on ranks](recipes/meshes.md#run-ranks-mrf) | [Meshes](recipes/meshes.md) |
| <a id="structure-contact"></a>302 | [Let a reed touch its seat and its stop](recipes/structure.md#structure-contact) | [Structures](recipes/structure.md) |
| <a id="run-ranks-sliding"></a>304 | [Run a sliding rotor in conformal stepping on ranks](recipes/meshes.md#run-ranks-sliding) | [Meshes](recipes/meshes.md) |
| <a id="run-ranks-imprint"></a>305 | [Run a sliding rotor by imprinting on ranks](recipes/meshes.md#run-ranks-imprint) | [Meshes](recipes/meshes.md) |
| <a id="spray-evaporation"></a>306 | [Let the droplets evaporate](recipes/sprays.md#spray-evaporation) | [Lagrangian sprays](recipes/sprays.md) |
| <a id="structure-temperature-field"></a>308 | [Load a structure with the temperatures of a flow case](recipes/structure.md#structure-temperature-field) | [Structures](recipes/structure.md) |
| <a id="surface-check"></a>310 | [Check an STL surface before it becomes a body](recipes/immersed.md#surface-check) | [Immersed bodies](recipes/immersed.md) |
| <a id="spray-cooler-chat"></a>316 | [Ask an assistant an engineering question: size a spray cooler](recipes/runs.md#spray-cooler-chat) | [Run control and health](recipes/runs.md) |

## Adding a feature

Every feature ships with a gallery example, case files and a recipe (D18,
as in PyAeroFlow); see [CONTRIBUTING.md](../CONTRIBUTING.md). Put the recipe on the topic page that fits, or start a new page when a new physics or topic begins. Number it with its pull request (a second recipe of the same pull request gets a letter: 54b), give it an `<a id>` anchor and a line in that page's table of contents, and name the case-file keys it uses. The two tables above are regenerated after the merge (`python -m tools.docs.indexes`); pull requests leave them alone. A new topic page also needs a row in the first table.
