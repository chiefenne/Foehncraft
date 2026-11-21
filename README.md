<p align="center">
  <img src="assets/logo/foehncraft.svg" alt="Foehncraft logo: a peregrine falcon" width="180" height="180">
</p>

<h1 align="center">Foehncraft</h1>

<p align="center">
  An open-source, industry-strength, general purpose CFD solver, written in Rust and built for modern GPUs.
</p>

---

> **Status: stage 6 done (RANS turbulence, heat transfer in forced
> convection, rotating machinery), release 0.6.** The stage
> plan and every decision are in [docs/requirements.md](docs/requirements.md);
> [What works today](#what-works-today) lists what is on `main` now. Stage 7 is
> under way: the energy equation of the fluid and steady conjugate heat
> transfer with solid regions, laminar and turbulent, are in; transient
> conjugate heat transfer comes next, buoyancy and variable density later.

[![CI](https://github.com/chiefenne/Foehncraft/actions/workflows/ci.yml/badge.svg)](https://github.com/chiefenne/Foehncraft/actions/workflows/ci.yml)
[![Release](https://img.shields.io/github/v/release/chiefenne/Foehncraft)](https://github.com/chiefenne/Foehncraft/releases)
[![Rust 1.89+](https://img.shields.io/badge/rust-1.89%2B-orange.svg)](Cargo.toml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

<table>
  <tr>
    <td colspan="3"><a href="docs/gallery/flow.md#cylinder"><img src="assets/readme/cylinder_vortex_street.jpg" alt="Vorticity of the Kármán vortex street behind a cylinder, Re = 100" width="900" height="167"></a></td>
  </tr>
  <tr>
    <td><a href="docs/gallery/flow.md#cavity"><img src="assets/readme/cavity_3d_polyhedra.jpg" alt="Streamlines of the 3D lid-driven cavity at Re = 1000 on 1.03 million polyhedra" width="300" height="200"></a></td>
    <td><a href="docs/gallery/meshes.md#tjunction"><img src="assets/readme/tjunction.jpg" alt="Streamlines of a jet in crossflow in a T-junction on 1.03 million tetrahedra" width="300" height="200"></a></td>
    <td><a href="docs/gallery/meshes.md"><img src="assets/readme/pin_fins.jpg" alt="Streamlines and wall shear stress in a staggered pin-fin array on 1.29 million tetrahedra" width="300" height="200"></a></td>
  </tr>
  <tr>
    <td><a href="docs/gallery/turbulence.md#impinging-jet"><img src="assets/readme/impinging_jet_temperature.jpg" alt="Temperature and streamlines of a turbulent round jet impinging on a heated plate" width="300" height="200"></a></td>
    <td><a href="docs/gallery/flow.md#turbulent-rotor"><img src="assets/readme/stirred_vessel_mrf.jpg" alt="Speed and streamlines in a turbulent stirred vessel with a five-blade rotor and three baffles, moving reference frame" width="300" height="200"></a></td>
    <td><a href="docs/gallery/flow.md#several-rotors"><img src="assets/readme/twin_rotors_vorticity.jpg" alt="Vorticity of two rotors with three and five blades turning at their own speeds on a sliding mesh" width="300" height="200"></a></td>
  </tr>
</table>

<sub>Real results of Foehncraft, each linked to its gallery entry: the Kármán vortex street (Schäfer &amp; Turek benchmark); the lid-driven cavity on a million polyhedra; a T-junction and a pin-fin array on Gmsh tetrahedra; a turbulent jet on a heated plate (k-ε-ζ-f₀); a stirred vessel in water with k-ω SST and the moving reference frame; two rotors on a sliding mesh.</sub>

## What works today

Everything below is on `main`, verified by the cases of the
[verification manual](docs/verification/manual.md) and shown in the
[gallery](docs/gallery.md).

- **Meshes.** Arbitrary polyhedral cells in 2D and 3D; the Gmsh v4, SU2, CGNS (HDF5),
  OpenFOAM polyMesh, industry `.msh` and VTK (`.vtk`, `.vtu`) readers behind one
  mesh-format interface (2D SU2 files extruded one layer), a VTU writer;
  a built-in test mesher (hexahedra, tetrahedra, prisms, polyhedra; graded
  segments, the cylinder O-grid, annuli for rotors); a quality report with
  automatic remediation of skewed cells (Gauss pressure gradient, limiter,
  bounded Rhie–Chow face velocity); partitioning with METIS and halos;
  VTKHDF output for ParaView.
- **Numerics.** Cell-centred finite volumes with weighted least-squares
  gradients, Venkatakrishnan(-Wang) limiters, linear upwind and bounded
  central convection, non-orthogonal correction, the time-step-independent
  Rhie–Chow interpolation; BDF2 and Crank–Nicolson; reproducible
  (gather-based) assembly in double precision.
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
  of non-circular ducts. Second order verified with manufactured
  solutions of every model on hex, polyhedral and tetrahedral cells
  (except ζ-f₀); validated on the flat plate, the channel against DNS, the backward-facing step and
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
  the same heat flow as the independent wall-function path to 10⁻⁴ (steady
  so far; transient conjugate heat transfer comes next).
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
  outlets.
- **Cases and properties.** TOML case files with checks and guidance;
  boundary values, sources and viscosities as expressions; a material
  database (air, water, oils, glycols) with sources.
- **Runs and results.** Restart, control files (`.stop`, `.result_N`,
  `.read_toml_N`), settings that change on a schedule, probes, forces,
  derived fields (vorticity, Q, λ₂, wall shear stress, y⁺, and the Reynolds stress, scales and wall values of the turbulence models),
  statistics of transient runs (time mean, RMS, covariance, phase averages,
  a confidence interval that says whether the window is long enough), mass and
  kinetic-energy ledgers (named terms per quantity, ready for further quantities), a divergence report, a log that can be
  circulated, JSON records for scripts and agents.
- **GPU and Python.** The flow loop on NVIDIA GPUs with results bitwise
  equal to the CPU; the Python package `foehncraft` that builds, runs and
  studies cases; the grid convergence index of every mesh sweep.

## Goal

Foehncraft aims to be a general purpose computational fluid dynamics solver of
industry strength: the class of tool that engineers rely on every day for
fluid flow, heat and mass transfer, species transport and reacting flows, on real
geometries and at real problem sizes. It is free and open source under the MIT
licence.

The guiding targets are:

- **General purpose.** Arbitrary polyhedral meshes in 3D, with a physics
  roadmap that reaches from incompressible flow to conjugate heat transfer,
  reacting flow, multiphase flow, moving meshes and scale-resolving
  turbulence.
- **Performance first.** Modern parallel GPUs are the primary platform, with
  a design target of more than 100 million cells on several GPUs and on
  clusters. The same solver runs on any CPU.
- **Accuracy that can be checked.** Double precision, numerics taken from the
  published literature, and every development stage verified against analytic
  solutions and published reference data.
- **Solver first.** A case file and a command line; meshes from established
  open tools, results to established open post-processors.

## First milestone

Version 1.0: incompressible flow, steady and transient, on unstructured 3D
meshes with the RANS models above, heat transfer with the energy equation,
buoyancy and conjugate heat transfer (stage 7), on NVIDIA GPUs and CPUs with
MPI (stage 8). Everything after that follows the roadmap in the
[requirements](docs/requirements.md): reacting flows, sprays, multiphase
flow, moving meshes and scale-resolving turbulence.

## Build and run

Needs Rust 1.89 or newer, a C compiler and CMake (HDF5 and METIS are built
from source on the first build; recipe 1 has the details).

```
cargo test --workspace
cargo build --release -p foehncraft-cli
target/release/foehn info
target/release/foehn bench g1 examples/cases/bench/g1_ci_small.toml
target/release/foehn mesh examples/cases/mesh/poly_partition.toml
target/release/foehn scalar examples/cases/scalar/mms_poly_box.toml
target/release/foehn poisson examples/cases/poisson/poisson_poly_box.toml
target/release/foehn flow examples/cases/flow/cavity_re100.toml
target/release/foehn flow examples/cases/flow/taylor_green_bdf2.toml
target/release/foehn flow examples/cases/flow/backstep_zeta_f0.toml
target/release/foehn flow examples/cases/flow/impinging_jet.toml
target/release/foehn flow examples/cases/flow/taylor_couette_mrf.toml
```

## Documentation

| Kind | What it is for | Where |
|---|---|---|
| Tutorial | a first mesh, end to end, in about ten minutes | [docs/tutorial.md](docs/tutorial.md) |
| How-to | steps for one task: build, run, check, partition, rerun G1, order of accuracy, pressure solver and hypre, GPU runs, flow cases and boundary conditions, transient runs, turbulence models and wall treatments, scalars and heat transfer, the coupled algorithm, rotating machinery, derived fields, restart and control files, scripts and agents, balances and divergence reports | [docs/recipes.md](docs/recipes.md) (index of topic pages) |
| Reference | every case file, its keys and defaults; the machine interfaces (JSON formats, exit codes, signals, diagnostic codes, balance and divergence files) | [docs/cases.md](docs/cases.md) (index of topic pages), [docs/keys.md](docs/keys.md), [docs/interfaces.md](docs/interfaces.md) |
| Examples | what each feature computes, with pictures and the commands | [docs/gallery.md](docs/gallery.md) (index of topic pages) |
| Engineering Gallery | industrial cases with complex geometry and the checks that trust them (demonstration entries, no reference data) | [docs/applications.md](docs/applications.md) |
| Sayings | the small words of the trade, for example *Residuen-gaffen* | [docs/sayings.md](docs/sayings.md) |
| Contributing | the documentation every feature ships with (gallery example, case files, recipe) and where it goes | [CONTRIBUTING.md](CONTRIBUTING.md) |
| Explanation | requirements, decisions, benchmark and verification results and their reading | [docs/requirements.md](docs/requirements.md), [G1](docs/benchmarks/g1.md), [stage 2](docs/verification/stage2.md), [stage 3](docs/verification/stage3.md), [stage 4](docs/verification/stage4.md), [stage 5](docs/verification/stage5.md), [stage 6](docs/verification/stage6.md), [verification manual](docs/verification/manual.md) |

The [changelog](CHANGELOG.md) lists every release.

| Crate | Content |
|---|---|
| `foehncraft-core` | generic floating-point type, index types, reproducible sums, counter-based random numbers (`no_std`) |
| `foehncraft-kernels` | compute kernels, one source for CPU and GPU (`no_std`): benchmark, scalar operators, multigrid, flow |
| `foehncraft-parallel` | communicator (serial, MPI), reproducible threaded reductions |
| `foehncraft-mesh` | polyhedral mesh, geometry, quality report, mesh-format layer with the Gmsh, SU2, CGNS, OpenFOAM, industry .msh and VTK readers, test mesher (also the cylinder benchmark mesh and annuli for rotors), partitioning and halo, VTKHDF output, point location, wall distance, the imprint of sliding interfaces |
| `foehncraft-fv` | finite-volume operators on a scalar (gradients, limiters, convection, diffusion, assembly), BiCGStab, and the pressure solver: CG with smoothed-aggregation multigrid in mixed precision, optional hypre reference; CPU backend (the GPU backend is in `foehncraft-cuda`) |
| `foehncraft-flow` | incompressible flow: SIMPLEC and SIMPLE (steady), PISO with BDF2 and Crank–Nicolson (transient), the coupled algorithm, the Rhie–Chow interpolation, boundary conditions and their check, turbulence models and wall treatments, passive scalars and thermal wall functions, the moving reference frame and sliding meshes, torque and machine performance, restart and control files, probes and force coefficients, the flow device written once over a CPU or GPU executor, the `foehn flow` case, exact solutions, reference and benchmark data |
| `foehncraft-bench` | gate G1: kernel benchmark, CPU backends |
| `foehncraft-cuda` | all NVIDIA GPU code (D60): the flow executor, the scalar and multigrid backend, the G1 backends with the CUDA C reference |
| `foehncraft-cli` | the `foehn` command |
| `foehncraft-py` | the Python package `foehncraft` (in `python/`) |
| `gpu/ptx-kernels` | the kernels compiled to NVIDIA PTX (nightly toolchain) |

## How it is built

Foehncraft follows the paradigm of the
[ArtCraft Crafting Apps](https://getartcraft.com/apps/craft-anything): open
source software in Rust, one engine with thin front ends on top, and
development carried out together with a coding agent in a loop of exploring,
building, checking and refining. The falcon joins the animals of that family.

What the agent builds is decided by a person. The development is driven by the
expert knowledge of the maintainer, from 40 years in computational fluid
dynamics: the requirements, the choice of methods, the test cases and the
judgement of every result come from that experience. The agent writes and
checks code; it does not decide what a CFD solver has to be.

Three rules follow from this:

1. **Requirements before code.** Each decision is recorded with its reason in
   [docs/requirements.md](docs/requirements.md) before anything is implemented.
2. **Every stage ends with results that can be tested.** No stage is done
   until its verification cases pass.
3. **Clean-room numerics.** Methods come from the published literature, with
   references; no code is taken from other solvers.

## Why a falcon

The peregrine falcon is the fastest animal there is, and its shape is the work
of the flow around it.

## Why the name

**Foehn** (German *Föhn*) is pronounced [føːn], with a long "ö" as in French *feu*, followed by "n". English speakers say /fɜːn/ ("fern") or /feɪn/. Foehncraft is therefore [ˈføːnkraft].

The foehn is the warm, dry wind that falls down the lee side of a mountain
range. Innsbruck, where the maintainer grew up and where the name comes from,
is one of the best-known foehn places in the world. With a southerly flow, air crosses the main Alpine ridge
at the Brenner Pass, the lowest gap at about 1370 m, runs down the Wipp valley
past Mount Patscherkofel and breaks into the Inn valley at Innsbruck, often in
strong gusts. Sinking without cloud, the air warms by about 1 °C per 100 m of
descent; when the foehn reaches the valley floor at about 580 m, the
temperature in town can rise by 5 to 10 °C within an hour, and the cold air
pooled in the valley is swept away.

Gap flow, a downslope windstorm, stratified flow over complex terrain and a
moving front between warm and cold air: the foehn is computational fluid
dynamics in the open air, and it has been studied at the University of
Innsbruck since the early twentieth century. The command is `foehn`.

## Contributing

The project is at stage 6 and has one maintainer. Comments on the
requirements are welcome as GitHub issues. [CONTRIBUTING.md](CONTRIBUTING.md)
has the rules every pull request follows, among them that it keeps this
README and the documentation up to date with what it changes.

## Licence

Distributed under the MIT licence. See [LICENSE](LICENSE) for more information.

<!-- year:start -->2026<!-- year:end --> Andreas Ennemoser – andreas.ennemoser@aon.at
