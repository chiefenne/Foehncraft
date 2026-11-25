# Solids moved by the fluid: spring valves and reed valves (plan, L52 and L53)

Plan, 2026-11-25. Requested by the maintainer in the immersed boundary
thread: "We leave in and out for the time being. What we add to the reqs is
a solid moving by fluid pressure. Eg a valve and fluid pressure supported by
a spring. The latter and damping modeled through a solver extension (no 3d
model) that could be entered by the user as well. Also reed valves shall be
able to be calculated. That means some small linear structural solver (but
make it general enough) to be able to run reed valves under fluid pressure.
The motion of reed valves also depends how stiff they are mounted. Research
all and let's do our card game."

Documentation only; **coding of each stage starts on the maintainer's
explicit go.** Bodies that enter or leave the domain stay out of scope (the
maintainer's "in and out" remark). State read: `main` at 5605686,
`docs/requirements.md` (D23, D30, D31, D53, D56, D58, D61, D91, D112, D119,
D122 (3), D182, D186, D201, D226, D234, L4, L13, L19, L22, L36), the
[L13 plan](../requirements/l13-immersed-boundary-plan.md) (stage 5, D234)
and the crates `foehncraft-flow` (`unsteady.rs`, `transient.rs`,
`sliding.rs`, `immersed.rs`, `restart.rs`, `balance.rs`, `monitor.rs`),
`foehncraft-mesh` (`geometry.rs`, `imprint.rs`) and `foehncraft-core`
(`predicates.rs`).

Clean room (D201): methods from papers and textbooks; measurements, analytic
solutions and correlations as references. No copyleft source is read (this
covers the open structural codes under the GPL as well as the CFD codes).
Numbers computed by other codes, including the reference values of the FSI
benchmarks, need a decision of their own (card 7). Commercial products are
named only where a feature is the published industry practice, never as a
source of code.

New ids: requirements **L52** (solids moved by fluid forces through a lumped
model) and **L53** (small linear structural solver and two-way FSI for reed
valves) are entered in `requirements.md` with this page. Decision numbers
**D286 to D292** are reserved for the cards of section 9 (provisional:
re-checked against `main` before the answers are recorded). Stage ids use
the prefix **MS** (moving solids).

## 1. Answer in short

1. **Everything rests on a deforming mesh, which is not built yet.** D91 has
   been a design constraint since stage 5, but the transient code still
   copies the cell volumes once at start-up, and there is no swept-volume
   face flux for deforming cells (only the rigid turn of the sliding rotors
   of D182, which keeps volumes). Stage **MS0** builds it and delivers the
   L22 piston (the evenly compressed cylinder) on the way, so that request
   from 2026-11-15 is answered by the same work.
2. **A valve on a spring is a rigid body with one degree of freedom plus a
   user-defined "mechanical model".** The fluid gives the force (pressure and
   shear integrated over the body's faces); the model gives spring, preload,
   damping, friction, stops and any extra force, from built-in elements,
   expressions (D53), tables and a linear state-space block. No 3D structure.
   This is **L52**, stages MS0 to MS4.
3. **A reed valve needs an elastic structure**: a small linear finite-element
   solver with beam, shell and solid elements, modal analysis, mounting
   springs, contact with the seat and the stop, and generalised-α time
   stepping, in its own crate. It is general enough for flaps, membranes,
   check-valve discs and later the thermal stress of L36. This is **L53**,
   stages MS5 and MS6.
4. **The coupling is the hard part, not the solvers.** A light reed in a
   dense fluid makes the plain staggered coupling unstable at any time step
   (added-mass effect). We plan strong coupling inside the time step with
   interface quasi-Newton acceleration (IQN-ILS) by default, Aitken as the
   simple fallback, and for the one-degree-of-freedom body an added-mass
   estimate that makes the scalar coupling stable even for light bodies.
5. **Valve closure** (the gap goes to zero) is handled by a **minimum gap**
   kept by the mesh motion and by **sealing faces** that switch the gap
   faces to walls while the valve is closed. Contact itself lives in the
   mechanical or structural model (stops with restitution, unilateral
   contact). Brinkman penalisation (#270) is the second option for the
   sealed gap, but the coupled solver breaks down at large μ/K (memory of
   #270), so face switching is the default.
6. **Sizing** (agent time, CI separate): MS0 to MS4 (spring valve) code
   about 6 to 9 h, runs 4 to 7 h; MS5 and MS6 (reed valve) code about 6 to
   9 h, runs 4 to 8 h; MS7 (free 6-DOF) and the compressible choked valve
   are **not sized yet**.

## 2. What the two requests ask for

### 2.1 Solid moved by fluid pressure, spring and damper (L52)

Typical devices: spring-loaded check valves and relief valves (poppet or
disc), pressure regulators, swing check valves (a flap on a hinge: one
rotational degree of freedom), flow-control pistons, ball check valves
(free body in a cage, later). What matters for the user:

- lift history x(t), opening and closing times, chatter and flutter
  (self-excited oscillation, the classic relief valve problem), the flow
  rate against the pressure difference, the force on the seat;
- the mechanical model entered by the user, not meshed: spring rate and
  preload, progressive springs (tables), viscous and Coulomb damping, the
  moving mass including the spring's effective mass (a third of the spring
  mass for a uniform coil spring), stops with restitution, an actuator force
  (solenoid) as an expression of time or of a monitored value;
- the force on the body from the fluid: pressure and viscous stress on all
  of its faces, including the faces in the gap.

### 2.2 Reed valves (L53)

Typical devices: compressor suction and discharge reeds (thin spring steel,
0.1 to 0.5 mm, lift a few mm), two-stroke and pump inlet reeds (steel or
composite petals), flapper valves, flexible flaps and membranes. What
matters:

- the reed is a thin plate: bending dominates, the first two or three modes
  shape the motion, the lift is large compared with the thickness but the
  rotations stay small (linear theory is the standard model in compressor
  valve design, Soedel 1992);
- **the mounting stiffness changes the motion** (maintainer): a clamp is
  never perfectly rigid, and the free length of the reed changes when it
  rolls off a curved stop or a clamp edge. We model it as distributed
  translational and rotational springs along the clamp line, with the clamp
  parts optionally meshed as solid elements in contact (section 6.4);
- contact with the seat (closing, impact, sticking by an oil film later)
  and with the stop or backer (limits lift, gives a second stiffness);
- the gas flow through the opening gap is fast and in compressors often
  near choked: this needs the pressure-based all-speed form noted in L22,
  which does not exist yet (section 8, risk R1). Liquid and low-speed gas
  reeds do not need it.

## 3. What exists on `main` (5605686)

| Piece | State | Used for |
|---|---|---|
| Rigid zone motion and grid flux | sliding rotors D182, D186: the grid flux ρ w·S of moving faces, imprinted interfaces per step, wall distance per step | the pattern for moving faces; a translating sleeve is a surface of revolution, so the D186 imprint already covers a poppet sliding in a bore |
| Deforming meshes (D91) | design constraint only: the ledger takes volumes per evaluation and the restart file has optional `points_*` and `volume_<k>` datasets, but the transient code copies the volumes once and no face carries a swept volume | MS0 builds it |
| Transient solver | PISO with outer correctors, BDF1, BDF2 with variable steps, Crank–Nicolson; Courant number check | the time loop that the coupling iterates in |
| Forces on boundaries | monitors D43 (force coefficients), immersed force Σ s u (#270) | the fluid force on a body |
| Penalisation | `[immersed]` Brinkman sink (#270) | option for the sealed gap and for a body in a cage |
| Exact predicates | `foehncraft-core::predicates` (orient2d, orient3d) | contact and inversion checks |
| Expressions, tables (D53) | strings and Python objects, same form | user laws for the mechanical model |
| Ledger as data (D119) | modules book rows per account | a mechanical-energy account: work of the fluid on the body = kinetic + spring energy + dissipation |
| Solid regions (D112), CHT | solid cells and interface faces | a heated valve later; the structural solver of L53 is separate (thin reeds do not fit a CFD solid region) |
| Not there | mesh motion solver, rigid-body dynamics, coupling iterations, any structural solver, contact, the all-speed pressure correction (L10/L22), cavitation (L3) | – |

## 4. How established codes and papers do it (research)

### 4.1 Coupling: partitioned, strong, accelerated

- **Partitioned coupling is the industry practice** for valve FSI: a flow
  solver and a structure solver exchange interface forces and displacements
  (commercial tool chains do it through a coupling service between their
  flow and structure products; the open coupling library preCICE does it
  between independent codes, Bungartz et al. 2016, Chourdakis et al. 2022).
  Monolithic solvers exist in research (Hron & Turek 2006) but are rarely
  used for valves.
- **Added-mass instability.** For an incompressible fluid the staggered
  (explicit) scheme is unstable whenever the structure is light compared
  with the fluid it moves, and a smaller time step makes it worse, not
  better (Causin, Gerbeau & Nobile 2005; Förster, Wall & Ramm 2007). A
  steel reed in air is mostly safe; a disc or reed in water or oil is not.
- **Strong coupling with acceleration.** Sub-iterations inside the time step
  until the interface converges, accelerated by dynamic Aitken relaxation
  (Küttler & Wall 2008) or by interface quasi-Newton methods: IQN-ILS
  (Degroote, Bathe & Vierendeels 2009), and multi-vector forms with linear
  cost (Lindner et al. 2015; Spenke, Hosters & Behr 2020). IQN-ILS is the
  default in the open coupling library and in most published valve work.
- **Added-mass partitioned (AMP) schemes** (Banks, Henshaw & Schwendeman
  2014 for beams; Banks, Henshaw, Schwendeman & Tang 2017 and 2018 for rigid
  bodies in incompressible flow; Banks, Henshaw & Sjögreen 2013 for light
  rigid bodies in compressible flow): an added-mass tensor, computed from
  the pressure response to a unit acceleration of the body, enters the
  body's equation and the pressure boundary condition, so the staggered
  scheme stays stable even for bodies of zero mass, without sub-iterations.
  For a body with one degree of freedom this is one extra pressure solve per
  mesh state (a scalar). We use it as a stabiliser of the strong coupling:
  the scalar added mass m_a goes on the left of the body equation and
  m_a·a^(k) on the right, which keeps the iteration contracting for light
  bodies and converges to the same solution.

### 4.2 Mesh motion around valves

- **Morphing (deforming mesh)** for strokes of a few cell sizes: Laplace
  smoothing with variable diffusivity (Jasak & Tuković 2006), the
  pseudo-elastic analogy with Jacobian-based stiffening (Stein, Tezduyar &
  Benney 2003), radial basis function interpolation of the boundary motion
  (de Boer, van der Schoot & Bijl 2007; greedy point selection, Rendall &
  Allen 2008), inverse distance weighting (Witteveen & Bijl 2009).
- **Layering and remeshing** for long strokes (engine valves, pistons): cell
  layers added or removed at a moving face; local remeshing when the quality
  falls (commercial dynamic-mesh tools offer smoothing, layering and
  remeshing together). Our O16 discussion was meant for this.
- **Rigid zone plus a deforming zone plus a sliding interface**: the poppet
  and its neighbourhood move rigidly inside a bore, only the thin region in
  front of the seat deforms. We have the sliding imprint for surfaces of
  revolution (D186), so this is close.
- **Overset meshes**: common in commercial tools for valves with large
  motion, and some offer a "zero gap" treatment that blocks the flow where
  two walls meet. Overset is not on our roadmap except as L4 "later".
- **Immersed bodies** (our L13): the cut cells of stage 2 and the moving
  bodies of stage 5 (D234) could carry a poppet, but thin reeds are refused
  by Q4 of L13, so reeds need a body-fitted moving mesh.

### 4.3 Valve closure and contact

Published valve FSI work keeps a **minimum gap** (a few µm to tens of µm)
and blocks or nearly blocks the flow there, by a porous or blocked zone or
by switching the gap cells or faces off, because a gap of zero changes the
mesh topology. Commercial tools describe "contact detection" that turns the
gap cells into a flow barrier, and "zero gap" interfaces for overset meshes.
The structure side handles contact by penalty or Lagrange multipliers with a
rigid seat and stop; impact is modelled by a coefficient of restitution in
lumped models.

### 4.4 Reed valve models (compressor literature)

The open proceedings of the International Compressor Engineering Conference
(Purdue, docs.lib.purdue.edu/icec) are the main open source. Classic
practice (Soedel 1992; Matos, Prata & Deschamps 2002 and later papers of the
Deschamps group): the reed as a plate solved by finite elements and reduced
to a few modes; the flow either by effective force and flow areas (lumped)
or by CFD on a moving mesh (radial diffuser with a flexible reed). Recent
papers couple a commercial flow solver with a commercial or own FE solver,
strongly coupled. Mounting stiffness and the stop are known to change the
first natural frequency and the closing delay. Six of these papers were read
for this plan; the table in section 12 lists what each one does.

### 4.5 Relief and check valves with springs

Direct spring-operated safety relief valves with CFD and a one-degree-of-
freedom body on a moving mesh, validated against measured lift and pressure:
Beune, Kuerten & van Heumen 2012 (high-pressure gas, opening), Song et al.
2014 (relief valve on a pressure vessel). The stability of such valves
(chatter) is analysed with lumped models and quarter-wave pipe dynamics
(Hős, Champneys, Paul & McNeely 2014).

### 4.6 Computer-graphics candidates (owner rule: always benchmark them)

| Need | Classical method | Graphics candidates to benchmark |
|---|---|---|
| Mesh motion from boundary motion | Laplace with diffusivity, pseudo-elastic, RBF, IDW | as-rigid-as-possible deformation (Sorkine & Alexa 2007) on the volume mesh; harmonic and bounded biharmonic weights (Joshi et al. 2007; Jacobson et al. 2011); regularized Kelvinlets, closed-form elastic displacement fields (de Goes & James 2017); grid-free harmonic extension by walk on stars (Sawhney, Miller, Gkioulekas & Crane 2023; Robin conditions: Miller et al. 2024) |
| Keeping cells valid near closure | quality limits and remeshing | barrier functions that keep every element positive (locally injective maps, Schüller et al. 2013); foldover-free untangling (Garanzha et al. 2021) as a repair step |
| Contact of the reed with seat and stop | penalty, Lagrange multipliers | Incremental Potential Contact (Li et al. 2020), intersection-free by a log barrier with continuous collision detection, and its codimensional form for thin shells with thickness (Li, Kaufman & Jiang 2021); exact continuous collision detection (Wang et al. 2021 benchmark) |
| Gap and self-contact of surfaces | proximity search | repulsive energies (Repulsive Shells, Sassen et al. 2024; Untangling Surfaces, Minarčík, Liu, Crane & Li 2026) to repair a surface that has crossed |

Read for this plan: the abstracts on cs.cmu.edu/~kmcrane (walk on stars,
walkin' Robin, repulsive shells, untangling surfaces, differential walk on
spheres). The full papers are read when the stage starts. Expectation
(to be measured, not assumed): RBF and the stiffened Laplace win on accuracy
per cost for small strokes; ARAP and biharmonic weights give better cell
shapes for large rotations of a flap; Kelvinlets are attractive as a cheap
closed-form predictor on GPU; walk on stars is noisy point by point, so its
cells need many walks before they are smooth enough (likely loses for mesh
motion, but it needs no volume solve); IPC is the strongest candidate for
reed contact because it never lets the reed pass the seat and gives a
smooth force for the coupling iteration.

## 5. L52: the lumped moving body (stages MS0 to MS4)

### 5.1 Body and mechanical model

A **moving body** is a set of wall patches (or an immersed body later) that
moves rigidly with q(t), one degree of freedom first: translation along an
axis or rotation about a hinge axis. The equation is

M q̈ + F_model(q, q̇, t, s) = F_fluid(q, q̇, t) + F_ext(t)

with F_fluid the projection of ∫(−p n + τ·n) dA over the body's faces on the
degree of freedom (a force or a torque) and s the internal states of the
model. The **mechanical model** is a small block diagram the user enters in
the case file, with these parts:

1. **built-in elements**: mass (with an effective spring mass), linear
   spring with preload, spring from a table (progressive, two-stage),
   viscous damper, regularised Coulomb friction (tanh law, so the step stays
   smooth), stops (lower and upper, penalty or exact with a coefficient of
   restitution), seat with sealing (section 5.3);
2. **expressions** (D53) of q, q̇, t, monitored values (a probe pressure, a
   flow rate) and states, for any force law the elements do not cover;
3. **tables** (D53) of force against lift;
4. a **linear state-space block** (A, B, C, D matrices) for small linear
   systems: a spring with surge modes, a two-mass pilot valve, a pipe
   resonator in lumped form;
5. later options, decided by card 3: an FMI 3.0 co-simulation unit (the open
   standard for system models from Modelica and similar tools) and
   sandboxed user code.

Time integration of the body: generalised-α (Chung & Hulbert 1993) with
the same Δt as the flow; for the linear block the exact matrix exponential.
The body state is part of the restart file (D42, W6), written on rank 0 and
broadcast, so runs on ranks stay bitwise (D31).

**Outputs** per step: q, q̇, q̈, the fluid force split into pressure and
shear, model forces per element, contact state, coupling iterations and
residual; in `run.jsonl`, `result.json` and a `body_<name>.csv`. A
**mechanical-energy account** in the ledger (D119): work of the fluid on the
body = Δ(kinetic + spring energy) + dissipated energy (damper, friction,
impacts), closed per step to the coupling tolerance.

### 5.2 Mesh motion (MS0, MS1)

- **MS0 deforming mesh core** (D91 made real): points and geometry per
  step on host and device, cell volumes per time level, swept-volume face
  fluxes from the face's swept volume so that the space conservation law
  holds exactly (Demirdžić & Perić 1988), BDF1 and BDF2 with variable steps
  on moving cells, PISO and the coupled solver, restart with points and old
  volumes, ranks (point halos) and GPU. First gallery case: the L22 evenly
  compressed cylinder (prescribed motion), with the energy equation and
  variable density (7c) for the adiabatic compression check.
- **MS1 mesh motion solver**: the displacement of the moving walls
  propagates into the mesh by one of the methods of section 4.2, chosen by
  the benchmark (card 1). Zones: rigid (moves with the body), deforming
  (morphing), fixed. Quality guard every step: minimum scaled Jacobian and
  face orthogonality with the exact predicates, the run stops with a clear
  message (a D30-style code) before a cell inverts, and names the place.

### 5.3 Closure (MS3)

- The mesh motion keeps the gap at g_min (case key, default a small
  fraction of the stroke); the body's lower stop sits at g_min.
- **Sealing faces**: the user names the seat line (a seat patch pair or a
  seat ring); when the body reaches the stop, the internal faces of the gap
  that cross the seat switch to walls (zero mass flux, exact), and switch
  back when the model's net force opens the valve. Pressure on the body's
  seat area is taken from the sealed side, so the opening force is right.
- Option: Brinkman blocking of the gap cells (#270) for comparison.
- The time step is cut around impacts (an event detector on the stop
  distance), and the coupling restarts its quasi-Newton history at each
  contact change.

### 5.4 Coupling (MS2)

Inside every time step: (1) predict q from the body equation, (2) move the
mesh, (3) solve the flow (outer correctors), (4) integrate F_fluid,
(5) update q from the body equation with the added-mass stabiliser,
(6) repeat until |Δq| and the force change are below tolerance (default
1e-6 of the stroke) or the iteration limit; acceleration by Aitken for one
degree of freedom, IQN-ILS for several. The explicit (one pass) scheme is
allowed only when the mass ratio check says it is stable (it warns, MS-code
in D30 style, and refuses for a liquid with a light body).

### 5.5 Verification and validation (MS4)

| Case | Reference | Check |
|---|---|---|
| Free-stream preservation on a deforming mesh (MS0) | exact | uniform flow stays uniform to round-off |
| Piston, slow adiabatic compression (MS0, L22) | p V^γ constant | entropy rise → 0 with Δt and h |
| Spring-mass in vacuum (no fluid) | analytic | frequency and energy to round-off |
| Cylinder or piston accelerating in still fluid | analytic added mass (½ρV for a cylinder in unbounded fluid; piston in a closed tube) | force history, light-body stability |
| Elastically mounted cylinder at low Re (vortex-induced vibration) | published lock-in range and amplitude at Re = 100 (laminar, 2D; sources checked when the stage starts) | amplitude and frequency against reduced velocity |
| Spring-loaded relief valve opening (gas or water) | measured lift and pressure (source with a measured lift history to be chosen when the stage starts; Song et al. 2014 measured blowdown, Beune et al. 2012 gas valves at 40 bar) | lift history, opening time, flow |
| Check valve closure and chatter | lumped stability limits (Hős et al. 2014) | onset of chatter |

Gallery: a spring-loaded disc valve opening and closing, speed and pressure
fields at three lifts (no four-fold or symbol-like patterns: an
axisymmetric wedge or a single offset port).

## 6. L53: the linear structural solver and reed valve FSI (stages MS5, MS6)

### 6.1 Scope of the solver

A new crate `foehncraft-structure`, small and general:

- **linear elastodynamics, small strains**, isotropic and orthotropic
  materials (the L48 material frames reused);
- **elements**: Euler–Bernoulli and Timoshenko beams (2D reed models, quick
  checks), Reissner–Mindlin shells with assumed strains against locking
  (MITC4, Dvorkin & Bathe 1984; MITC3+ for triangles, Lee & Bathe 2014) for
  reeds and flaps, quadratic tetrahedra and hexahedra for clamps, thick
  plates and general parts;
- **analyses**: static, modal (shift-invert Lanczos on K − ω²M, the first
  modes for checks and for reduction), transient by generalised-α, mode
  superposition (reduced model of the first n modes) and Craig–Bampton
  reduction (1968) of a clamp substructure;
- **damping**: Rayleigh and modal damping ratios;
- **boundary conditions**: fixed, prescribed motion, point and distributed
  springs (translational and rotational per unit length or area: the
  mounting stiffness), pressure and traction loads, follower loads from the
  fluid;
- **contact**: unilateral contact of nodes against rigid seat and stop
  surfaces (penalty first, active set with exact gap second; IPC as the
  graphics candidate, section 4.6), frictionless first;
- **meshes**: shell meshes from Gmsh (`.msh`) or built from a reed outline
  (length, width, thickness, port holes) with the quad mesher in Gmsh;
  solids from Gmsh. The structural mesh is separate from the CFD mesh.
- **solvers**: the sparse direct factorisation (structural systems are
  small, up to about 10⁵ unknowns) with the existing iterative solvers as a
  fallback; CPU on rank 0, the flow on all ranks.

Geometric nonlinearity (large rotations, the flag of Turek & Hron FSI2 and
FSI3) is **not** in the first version; a corotational shell is the natural
next step and is noted as an open item.

### 6.2 Coupling with the flow

- **Interface mapping**: the fluid wall faces on both sides of the reed map
  to the shell mid-surface; displacements go from the shell to the fluid
  nodes by the shell's own shape functions plus the thickness offset, and
  forces go back by the transpose of the same operator, so the work done on
  the interface is identical on both sides (Farhat, Lesoinne & Le Tallec
  1998). Nearest-projection and RBF mapping are the alternatives (compared
  by de Boer, van Zuijlen & Bijl 2008).
- **Strong coupling** as in 5.4 with IQN-ILS on the interface displacement,
  history reused over a few time steps, restarted at contact changes.
- **Ledger**: the work on the interface (fluid side and structure side)
  books as an exchange pair; elastic and kinetic energy of the structure as
  storage; closes to the coupling tolerance.

### 6.3 Mounting stiffness

Three ways, from simple to complete:

1. **clamp springs**: translational and rotational springs per unit length
   along the clamp line (user values, or calibrated from a measured first
   natural frequency by a small fit tool);
2. **clamp parts as solids**: the clamp plates and the valve plate meshed
   with solid elements, the reed in frictionless contact with them, so the
   effective free length follows the load (the reed rolls off the clamp
   edge);
3. **stop or backer contact** as a second stiffness: the reed lies on a
   curved stop and its free length shortens as it lifts.

### 6.4 Verification and validation (MS6)

| Case | Reference | Check |
|---|---|---|
| Cantilever beam and plate, static and modal | analytic beam theory; cantilever plate frequencies (Leissa 1969, NASA SP-160) | deflection, first five frequencies |
| Shell patch tests and locking tests | MacNeal & Harder 1985 standard set; Scordelis–Lo roof (Belytschko et al. 1985) | convergence, no locking |
| Contact of a beam on a rigid stop | analytic | contact force, gap ≥ 0 |
| Clamp spring calibration | analytic beam with rotational end spring | frequency against clamp stiffness |
| Flexible flap or flag behind a cylinder, small deformation | Turek & Hron (2006) FSI1 and the structure tests CSM1 (small part) | tip displacement, drag and lift; needs D201 card 7 because these are computed numbers |
| Compressor reed valve with flow | measured reed lift and natural frequency (Gasche et al. as used by González et al. 2017: first mode 31.4 Hz, lift frequency 57.9 Hz; Lemke et al. 2016 flapper valve, 366 Hz) | lift history, opening and closing timing, natural frequency; target better than the 16 % peak-opening error of the published modal model |

Gallery: a reed opening over a port, deflection and speed fields at
several crank angles; a flap behind a cylinder (FSI1). Each with 1 to 3
colourful field images (working rule 8) and the mesh.

## 7. Interfaces (D77 to D87, D196)

Everything is case-file driven and scriptable, with structured outputs:

- `[[moving_body]]` with `patches`, `dof = "translate" | "rotate"`, axis,
  `[moving_body.model]` with element tables and expressions, `[mesh_motion]`
  with method and zones, `[moving_body.seal]` for the seat;
- `[[structure]]` with the shell or solid mesh, material, clamp springs,
  contact surfaces, and `[structure.coupling]` with the fluid patches;
- tools in the registry (D196): `structure.modes` (modal analysis alone),
  `structure.static`, `mesh_motion.check` (moves the mesh through the stroke
  without flow and reports the worst cell), `moving_body.estimate`
  (mass ratio, added mass, stable time step, explicit or strong coupling);
  each callable from the CLI, Python and MCP;
- the journal (D197) records the coupling iterations per step.

## 8. Risks

| Id | Risk | Effect | Mitigation |
|---|---|---|---|
| R1 | Gas valves with choked gap flow need the pressure-based all-speed form (L10, L22 note) | compressor reeds at real pressure ratios cannot be run | first cases liquid or low-speed gas; the all-speed form as its own stage after MS6 (card 6) |
| R2 | Added-mass instability for light bodies in liquid | divergence | strong coupling, added-mass stabiliser, refusal of the explicit scheme |
| R3 | Mesh quality near closure (thin gap cells, high aspect ratio) | stiff pressure solve, inverted cells | minimum gap, quality guard, sealing faces instead of penalisation, morphing benchmark |
| R4 | Impacts need very small time steps | long runs | event detection and step cut only around contact |
| R5 | Linear structure for large reed deflection | error in lift | check against the measured lift; corotational shell as the next step |
| R6 | Turbulence in the gap with moving walls | wall distance and wall functions per step | wall distance per step exists for sliding; the gap is often laminar or transitional (report y⁺ and Re of the gap) |
| R7 | Cavitation in liquid valves (L3 not built) | wrong forces at large Δp | first cases below cavitation; flagged in limitations (D199) |
| R8 | Coupled solver and moving meshes | the coupled preconditioner assumes fixed geometry in its setup | rebuild per step as for the imprinted sliding interface; measure the cost |
| R9 | Clean room for FE elements | open structural codes are GPL | elements from the papers only (MITC, generalised-α, Craig–Bampton), never read GPL source |

## 9. Decisions for the maintainer (cards, one at a time)

Each card is posted in the thread, the next only after the answer.

| Card | Id | Question | Recommendation |
|---|---|---|---|
| 1 | D286 | How should the mesh follow a moving valve? | deforming mesh (morphing) with a rigid zone and sliding imprint for long strokes; morphing method chosen by a benchmark of classical and graphics methods; overset and immersed later |
| 2 | D287 | How should a closing valve be handled? | minimum gap with sealing faces switched to walls; Brinkman blocking as the compared option |
| 3 | D288 | How does the user enter spring, damping and other forces? | built-in elements, expressions, tables and a linear state-space block first; FMI co-simulation later |
| 4 | D289 | Which structural solver for reeds? | own linear FE crate with beams, MITC shells and solids, modal reduction, contact; geometric nonlinearity later |
| 5 | D290 | How is the mounting stiffness modelled? | clamp springs first (with calibration from a natural frequency), clamp parts in contact second |
| 6 | D291 | Which fluids come first? | liquid and low-speed gas valves first; the all-speed form for choked gas valves as its own stage after MS6 |
| 7 | D292 | Which references count for validation (D201)? | measurements and analytic results; the Turek–Hron FSI1 values allowed as an accepted benchmark like Schäfer & Turek, by its own decision |

## 10. Stages and sizing

Agent time (owner rule): code in minutes per PR, runs in hours, CI
separate. Each stage ends with gallery, case and recipe entries (D18,
working rule 8).

| Stage | Content | PRs | Code | Runs |
|---|---|---|---|---|
| MS0 | Deforming mesh core (D91), L22 evenly compressed piston | 2–3 | 1.5–2.5 h | 1 h |
| MS1 | Mesh motion solver + zones + quality guard + morphing benchmark (section 4.6) | 2–3 | 1.5–2.5 h | 1–2 h |
| MS2 | Moving body with mechanical model, coupling (Aitken, added mass), restart, ledger | 2 | 1.5–2 h | 1 h |
| MS3 | Closure: minimum gap, sealing faces, impact events | 1–2 | 1–1.5 h | 1 h |
| MS4 | Validation: added mass, VIV, relief valve, gallery | 1–2 | 0.5–1 h | 1–2 h |
| MS5 | Structural crate: elements, modal, generalised-α, springs, contact, verification | 3–4 | 3–4.5 h | 1 h |
| MS6 | Reed FSI: mapping, IQN-ILS, ledger pair, flap and reed validation, gallery | 2–3 | 2–3 h | 3–6 h |
| MS7 | Free rigid body with up to 6 degrees of freedom (L4), ball valve | – | not sized yet | – |
| MS8 | All-speed pressure correction for choked valve flow (with L10) | – | not sized yet | – |
| MS9 | Corotational shells (large rotations), Turek–Hron FSI2/FSI3 | – | not sized yet | – |

Totals: L52 (MS0 to MS4) code 6 to 9 h, runs 4 to 7 h; L53 (MS5, MS6) code
5 to 7.5 h, runs 4 to 7 h. MS5 can run in parallel with MS1 to MS4 in its
own thread (no shared files except the case schema).

Order: MS0 → MS1 → MS2 → MS3 → MS4; MS5 in parallel from the start;
MS6 after MS3 and MS5. MPI ranks: MS0 must keep W7 working (geometry halos
per step); body and structure on rank 0 with broadcast.

## 11. Flaws of industry practice to avoid (for the challenger)

- **Silent explicit coupling.** Tools that let a user run a staggered
  coupling for a light body in water and blow up later; we estimate the
  mass ratio first and refuse or switch to strong coupling.
- **Porous blocking with a tuned resistance.** The sealed gap leaks by an
  amount that depends on a resistance and a cell size; our face switching
  is exact (zero flux) and reports the leak as zero by construction.
- **Interface maps that do not conserve work.** Interpolating forces and
  displacements separately creates energy at the interface; our transposed
  operator conserves it, and the ledger shows it.
- **Mesh motion that fails late.** Inverted cells found only when the
  solver diverges; our quality guard stops before inversion and
  `mesh_motion.check` runs the stroke without flow in seconds.
- **Clamp as perfectly rigid by default.** We ask for the mounting
  stiffness or use the calibration tool, and say in the result when the
  clamp was taken as rigid.

## 12. Literature read for this plan

Read on 2026-11-25 (full text unless marked). The numbers below are what the
papers report; a stage reads the paper again before using any of them.

| Paper | Structure | Mesh motion | Closure | Coupling | Validation |
|---|---|---|---|---|---|
| Estruch, Lehmkuhl, Rigola & Pérez-Segarra 2014, ICEC 2380 (reed under piston displacement) | plate by mode summation, up to 3 cantilever modes, Crank–Nicolson | RBF morphing (Wendland C2), grid velocity from the space conservation law | constraint w ≥ 0 and a restitution correction on the modal velocities; seat as an immersed body | loose (one pass per step), justified by ρ_f/ρ_s = 8.9e-4 | none against measurement |
| González, Lehmkuhl, Naseri, Rigola & Oliva 2016, ICEC 2490 | Kirchhoff–Love plate, 4 modes from a clamped FE modal analysis, Newmark | pseudo-elastic with stiffening of small cells; immersed boundary turns the seat region fluid as the gap opens | restitution 0.3 on the modal velocity | semi-implicit (pressure implicit), block Gauss–Seidel with Aitken | none; coarse mesh noted by the authors |
| González, Naseri, Rigola, Pérez-Segarra & Oliva 2017, IOP Conf. Ser. Mater. Sci. Eng. 232 012032 | as above | as above | as above | as above | against the measurements of Gasche et al. (Re = 10 000): first mode 31.1 Hz vs 31.4 Hz measured, lift frequency 58.6 vs 57.9 Hz, peak opening 16 % too high (blamed on the impact model) |
| Möhl, Thomas & Hesse 2016, ICEC 2390 | reed as an FE solid (commercial tool chain), 2D | sub-domains with stiffness by expressions | limiter at 2.5 mm, minimum gap 0.1 mm in the contact pair | iterative two-way, scheme not stated | discharge coefficient within 0.4 % of Böswirth (1994); lift "fair" against a 1D model |
| Lemke, König, Hennig, Försterling & Köhler 2016, ICEC 2415 (flapper valve) | FE solids (reed, retainer, seat) with contact, coupled with a commercial flow solver | overset mesh morphed with the valve, "zero gap" contact | contact in the FE model, zero-gap blocking in the flow | exchange every step, Δt = 1.5e-5 s | first mode 400 Hz vs 366 Hz measured (9 %); lift by laser vibrometer, frequencies agree, amplitudes low; the 1-DOF model misses the oscillations |
| Wang et al. 2021, Appl. Sci. 11, 3946 | 3D solid FE, contact with seat and limiter | not described | contact in the FE model | scheme not stated | cylinder pressure and suction lift (eddy-current sensor) "in good agreement", no metric |
| Song et al. 2014, Energy Convers. Manag. 81 (spring safety relief valve) | 1-DOF disc, m = 0.96 kg, k = 22.3 kN/m, no damping | displacement diffusion over four hex blocks | minimum lift 0.4 mm, no contact model | explicit (forces of the previous step), Δt = 1e-4 s | blowdown within 6 to 11 % of the measurement; lift history not measured |
| Beune, Kuerten & van Heumen 2012, Computers & Fluids 64 (abstract only) | 1-DOF by Newton's law | – | – | – | quasi-steady models miss the flow-history effect; oscillation on closing at 40 bar |
| Banks, Henshaw, Schwendeman & Tang 2016/2017 (arXiv 1611.05711) | rigid body | – | – | AMP: Robin pressure condition with the body acceleration, stable without sub-iterations for light and zero-mass bodies | second-order convergence, no measurement |
| Banks, Henshaw & Sjögreen 2012/2013 (arXiv 1208.3511) | rigid body with added-mass tensors, DIRK | overlapping grids | – | AMP for compressible flow, stable for zero mass | analytic receding piston |
| Spenke, Hosters & Behr 2020 (arXiv 2001.07947) | – | – | – | multi-vector quasi-Newton with linear cost | (test section not read) |

What this tells us:

- **The published reed models are modal plates with a restitution
  correction**, coupled loosely or semi-implicitly; the best documented
  agreement with measurement (González et al. 2017) gets the frequency to
  1 % and over-predicts the peak opening by 16 %, and the authors blame the
  impact model. A real FE plate with contact (MS5) and strong coupling
  (MS6) attacks exactly that weak point.
- **Every paper keeps a minimum gap** (0.1 mm, 0.4 mm) or blocks the gap
  (immersed seat, zero-gap overset). None reports the leak through the
  "closed" valve; our sealing faces make it zero by construction.
- **Explicit coupling is common** even for the relief valve (Song et al.
  2014); it works for gas because the added mass is small, and is the flaw
  section 11 names for liquids.
- **Measured data with lift histories** exist in the open compressor
  proceedings (Gasche et al. as used by González et al. 2017; Lemke et al.
  2016 for a flapper valve). The original measurement papers are read when
  MS6 starts; the relief valve case needs a source with a measured lift
  history (Song et al. 2014 measured blowdown only).
- Matos, Prata & Deschamps (2002) and Soedel (1992) are cited from the
  literature above as the classic design references; not read for this plan.

## 13. References

- Banks, J. W., Henshaw, W. D., Schwendeman, D. W. (2014). An analysis of a
  new stable partitioned algorithm for FSI problems. Part I: incompressible
  flow and elastic solids. *J. Comput. Phys.* 269.
- Banks, J. W., Henshaw, W. D., Schwendeman, D. W., Tang, Q. (2017). A
  stable partitioned FSI algorithm for rigid bodies and incompressible flow.
  Part I: model problem analysis. *J. Comput. Phys.* 343 (arXiv 1611.05711);
  Part II (2017); three dimensions (2018, arXiv 1711.04742).
- Banks, J. W., Henshaw, W. D., Sjögreen, B. (2013). A stable FSI algorithm
  for light rigid bodies in compressible flow. *J. Comput. Phys.* 245
  (arXiv 1208.3511).
- Belytschko, T. et al. (1985). Stress projection for membrane and shear
  locking in shell finite elements. *Comput. Methods Appl. Mech. Eng.* 51.
- Beune, A., Kuerten, J. G. M., van Heumen, M. P. C. (2012). CFD analysis
  with fluid–structure interaction of opening high-pressure safety valves.
  *Computers & Fluids* 64.
- Bungartz, H.-J. et al. (2016). preCICE – a fully parallel library for
  multi-physics surface coupling. *Computers & Fluids* 141; Chourdakis, G.
  et al. (2022). preCICE v2. *Open Research Europe* 2.
- Causin, P., Gerbeau, J.-F., Nobile, F. (2005). Added-mass effect in the
  design of partitioned algorithms for fluid–structure problems. *Comput.
  Methods Appl. Mech. Eng.* 194.
- Chung, J., Hulbert, G. M. (1993). A time integration algorithm for
  structural dynamics with improved numerical dissipation: the generalized-α
  method. *J. Appl. Mech.* 60.
- Craig, R. R., Bampton, M. C. C. (1968). Coupling of substructures for
  dynamic analyses. *AIAA J.* 6.
- de Boer, A., van der Schoot, M. S., Bijl, H. (2007). Mesh deformation
  based on radial basis function interpolation. *Computers & Structures* 85.
- de Boer, A., van Zuijlen, A. H., Bijl, H. (2008). Comparison of
  conservative and consistent approaches for the coupling of non-matching
  meshes. *Comput. Methods Appl. Mech. Eng.* 197.
- de Goes, F., James, D. L. (2017). Regularized Kelvinlets: sculpting brushes
  based on fundamental solutions of elasticity. *ACM Trans. Graph.* 36.
- Degroote, J., Bathe, K.-J., Vierendeels, J. (2009). Performance of a new
  partitioned procedure versus a monolithic procedure in fluid–structure
  interaction. *Computers & Structures* 87.
- Demirdžić, I., Perić, M. (1988). Space conservation law in finite volume
  calculations of fluid flow. *Int. J. Numer. Methods Fluids* 8.
- Dvorkin, E. N., Bathe, K.-J. (1984). A continuum mechanics based
  four-node shell element for general nonlinear analysis. *Eng. Comput.* 1.
- Estruch, O., Lehmkuhl, O., Rigola, J., Pérez-Segarra, C. D. (2014).
  Fluid–structure interaction of a reed type valve subjected to piston
  displacement. *Int. Compressor Eng. Conf.*, Purdue, paper 2380.
- Farhat, C., Lesoinne, M., Le Tallec, P. (1998). Load and motion transfer
  algorithms for fluid/structure interaction problems with non-matching
  discrete interfaces. *Comput. Methods Appl. Mech. Eng.* 157.
- Förster, C., Wall, W. A., Ramm, E. (2007). Artificial added mass
  instabilities in sequential staggered coupling of nonlinear structures and
  incompressible viscous flows. *Comput. Methods Appl. Mech. Eng.* 196.
- Garanzha, V. et al. (2021). Foldover-free maps in 50 lines of code. *ACM
  Trans. Graph.* 40.
- González, I., Lehmkuhl, O., Naseri, A., Rigola, J., Oliva, A. (2016).
  Fluid–structure interaction of a reed type valve. *Int. Compressor Eng.
  Conf.*, Purdue, paper 2490.
- González, I., Naseri, A., Rigola, J., Pérez-Segarra, C. D., Oliva, A.
  (2017). A fluid–structure interaction solver for the fluid flow through
  reed type valves. *IOP Conf. Ser.: Mater. Sci. Eng.* 232, 012032.
- Hős, C. J., Champneys, A. R., Paul, K., McNeely, M. (2014). Dynamic
  behavior of direct spring loaded pressure relief valves in gas service:
  model development, measurements and instability mechanisms. *J. Loss
  Prev. Process Ind.* 31.
- Hron, J., Turek, S. (2006). A monolithic FEM/multigrid solver for an ALE
  formulation of fluid–structure interaction with applications in
  biomechanics; Turek, S., Hron, J. (2006). Proposal for numerical
  benchmarking of fluid–structure interaction between an elastic object and
  laminar incompressible flow. In *Fluid–Structure Interaction*, LNCSE 53,
  Springer.
- Jacobson, A., Baran, I., Popović, J., Sorkine, O. (2011). Bounded
  biharmonic weights for real-time deformation. *ACM Trans. Graph.* 30.
- Jasak, H., Tuković, Ž. (2006). Automatic mesh motion for the unstructured
  finite volume method. *Trans. FAMENA* 30.
- Joshi, P. et al. (2007). Harmonic coordinates for character articulation.
  *ACM Trans. Graph.* 26.
- Küttler, U., Wall, W. A. (2008). Fixed-point fluid–structure interaction
  solvers with dynamic relaxation. *Comput. Mech.* 43.
- Lee, P.-S., Bathe, K.-J. (2014) and Lee, Y. et al. (2014). The MITC3+
  shell element. *Computers & Structures* 138.
- Leissa, A. W. (1969). *Vibration of Plates.* NASA SP-160.
- Lemke, N. C., König, M., Hennig, J., Försterling, S., Köhler, J. (2016).
  Transient experimental and 3D-FSI investigation of flapper valve dynamics
  for refrigerant compressors. *Int. Compressor Eng. Conf.*, Purdue, paper
  2415.
- Li, M. et al. (2020). Incremental potential contact: intersection- and
  inversion-free, large-deformation dynamics. *ACM Trans. Graph.* 39.
- Li, M., Kaufman, D. M., Jiang, C. (2021). Codimensional incremental
  potential contact. *ACM Trans. Graph.* 40 (arXiv 2012.04457).
- Lindner, F., Mehl, M., Scheufele, K., Uekermann, B. (2015). A comparison of
  various quasi-Newton schemes for partitioned fluid–structure interaction.
  *Proc. Coupled Problems.*
- MacNeal, R. H., Harder, R. L. (1985). A proposed standard set of problems
  to test finite element accuracy. *Finite Elem. Anal. Des.* 1.
- Matos, F. F. S., Prata, A. T., Deschamps, C. J. (2002). Numerical
  simulation of the dynamics of reed type valves. *Int. Compressor Eng.
  Conf.*, Purdue.
- Miller, B., Sawhney, R., Crane, K., Gkioulekas, I. (2024). Walkin' Robin:
  walk on stars with Robin boundary conditions. *ACM Trans. Graph.* 43.
- Minarčík, J., Liu, M., Crane, K., Li, M. (2026). Untangling surfaces via
  shape and mesh repulsion. *ACM Trans. Graph.* 45.
- Möhl, C., Thomas, C., Hesse, U. (2016). An investigation into the
  dynamics of self-acting compressor valves. *Int. Compressor Eng. Conf.*,
  Purdue, paper 2390.
- Rendall, T. C. S., Allen, C. B. (2008). Unified fluid–structure
  interpolation and mesh motion using radial basis functions. *Int. J.
  Numer. Methods Eng.* 74.
- Sassen, J. et al. (2024). Repulsive shells. *ACM Trans. Graph.* 43.
- Sawhney, R., Miller, B., Gkioulekas, I., Crane, K. (2023). Walk on stars: a
  grid-free Monte Carlo method for PDEs with Neumann boundary conditions.
  *ACM Trans. Graph.* 42.
- Schüller, C., Kavan, L., Panozzo, D., Sorkine-Hornung, O. (2013). Locally
  injective mappings. *Comput. Graph. Forum* 32.
- Soedel, W. (1992). *Design and Mechanics of Compressor Valves.* Purdue
  University (Ray W. Herrick Laboratories).
- Song, X. et al. (2014). A CFD analysis of the dynamics of a
  direct-operated safety relief valve mounted on a pressure vessel. *Energy
  Convers. Manag.* 81.
- Sorkine, O., Alexa, M. (2007). As-rigid-as-possible surface modeling.
  *Proc. Symp. Geometry Processing.*
- Spenke, T., Hosters, N., Behr, M. (2020). A multi-vector interface
  quasi-Newton method with linear complexity for partitioned fluid–structure
  interaction. *Comput. Methods Appl. Mech. Eng.* 361 (arXiv 2001.07947).
- Stein, K., Tezduyar, T., Benney, R. (2003). Mesh moving techniques for
  fluid–structure interactions with large displacements. *J. Appl. Mech.* 70.
- Wang, B., Ferguson, Z., Schneider, T., Jiang, X., Attene, M., Panozzo, D.
  (2021). A large-scale benchmark and an inclusion-based algorithm for
  continuous collision detection. *ACM Trans. Graph.* 40.
- Wang, Y., Wang, J., He, Z., Sun, J., Wang, T., Liu, C. (2021).
  Investigation on dynamic characteristics of the reed valve in compressors
  based on fluid–structure interaction method. *Appl. Sci.* 11, 3946.
- Witteveen, J. A. S., Bijl, H. (2009). Explicit mesh deformation using
  inverse distance weighting interpolation. *AIAA paper* 2009-3996.
