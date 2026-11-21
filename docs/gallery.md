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
| [Meshes](gallery/meshes.md) | the built-in test mesher in 2D and 3D, partitioning, the Gmsh reader and flow on poor Gmsh meshes | 3, 4, 5, 6, 26, 127, 129, 134, 138 |
| [Scalar transport and the pressure solver](gallery/scalar-pressure.md) | convection–diffusion order of accuracy, convection schemes and limiters, the pressure multigrid | 8, 9, 10, 11 |
| [Incompressible flow](gallery/flow.md) | steady and transient laminar verification and benchmarks: cavity, Poiseuille, Taylor–Green, cylinder, time-step control | 12, 13, 14, 15, 20, 78, 80, 81, 83, 85, 86, 87, 94, 95, 96, 98, 99, 100, 101, 102, 103 |
| [Passive scalars and the wall distance](gallery/scalars.md) | passive scalars carried by the flow, the bounded BDF2, the exact wall distance | 63 |
| [Energy equation and heat transfer](gallery/energy.md) | the energy equation of the fluid with temperature-dependent properties, viscous heating and thermal walls: a heated cavity, conduction with k(T), Couette flow heated by friction, conjugate heat transfer with solid regions, laminar and turbulent (a finned heat sink) | 154, 156, 158, 159, 160 |
| [Turbulence models](gallery/turbulence.md) | the RANS turbulence models of stage 6 with their wall treatments: k-ω SST on the TMR flat plate and in a channel | 68, 71, 73, 74, 88, 89, 90, 113, 115, 141, 140 |
| [Boundary conditions](gallery/boundary-conditions.md) | slip walls, axisymmetric wedges, inlet profiles from tables, backflow at pressure outlets, convective outlets, periodic modules | 21, 67 |
| [Expressions and materials](gallery/properties.md) | boundary values, body forces and viscosities as expressions, the material database, regions, walls with a mass flux | 23 |
| [Python](gallery/python.md) | cases built, run and studied from Python: expressions with Python operators, callbacks, study tables | 24 |
| [Verification](gallery/verification.md) | the discretisation uncertainty of mesh sweeps (grid convergence index) and the generated verification manual | 25 |
| [Post-processing](gallery/post-processing.md) | derived fields: vorticity, Q, λ₂, wall shear stress, y⁺, from the run or a restart file | 22, 108, 110 |
| [Run control, health and adjoint readiness](gallery/runs.md) | adjoint-ready kernels, runs that a program drives, balances and the divergence report, the log, the table columns, scheduled changes | 17, 18, 19, 27, 28, 29 |

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

## Adding a feature

Every feature ships with a gallery example, case files and a recipe (D18,
as in PyAeroFlow); see [CONTRIBUTING.md](../CONTRIBUTING.md). Put the example on the topic page that fits, or start a new page when a new physics or topic begins (energy, turbulence, reacting flows, sprays, ...). Number it with its pull request (a second example of the same pull request gets a letter: 54b), give it an `<a id>` anchor and a row in that page's table of contents. The two tables above are regenerated after the merge (`python -m tools.docs.indexes`); pull requests leave them alone. A new topic page also needs a row in the first table.
