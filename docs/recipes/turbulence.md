# Turbulence models

Part of the [Recipes](../recipes.md): run a case with a RANS turbulence model, choose the wall treatment, set inlet turbulence and roughness, and check the near-wall mesh. Every key is on the [key page](../cases/flow-keys.md#turbulence-d129d139).

| # | Recipe |
|---|---|
| 68 | [Run a case with the k-ω SST model](#sst) |
| 71 | [Run a case with the k-ε family and wall functions](#k-epsilon) |
| 73 | [Run a case with k-ε-ζ-f and resolved walls](#zeta-f) |
| 74 | [Run a case with k-ε-ζ-f₀ at any first-cell y⁺](#zeta-f0) |
| 88 | [Run a case with the elliptic-blending k-ε](#elliptic-blending) |
| 89 | [Validate the turbulence models on the backward-facing step](#step-validation) |
| 90 | [Feed a jet with a developed pipe profile and compute the Nusselt number](#impinging-jet) |
| 113 | [Make SST sensitive to curvature and system rotation](#curvature-correction) |
| 115 | [Get the secondary flow of non-circular ducts with SST](#qcr) |
| 141 | [Choose the sink of k of the ζ-f models](#k-sink) |
| 140 | [Say what turbulence enters through an outlet](#outlet-backflow) |
| 194 | [Run the QCR on tetrahedra and perturbed meshes](#qcr-consistency) |
| 197 | [Run the curvature correction on tetrahedra and perturbed meshes](#curvature-consistency) |

---

<a id="sst"></a>
## 68. Run a case with the k-ω SST model

Switch the model on and pick the wall treatment by the mesh: `low-re`
when the first cell lies at y⁺ ≲ 1, `automatic` (the default) for any
other mesh, with the first cell anywhere from the sublayer to the log
layer (C_f of the flat plate within 2.6 % for y⁺ 0.24 to 85; TURB001
warns above y⁺ 300):

```toml
[turbulence]
model = "k-omega-sst"
wall_treatment = "low-re"        # or "automatic"
# production_limiter = "menter"  # "kato-launder", "durbin", "none"
# ambient = { k = 2.25e-7, omega = 125.0 }   # sustaining terms, external flows
```

Give the inlets their turbulence: an intensity with a length scale or a
viscosity ratio, or k and ω directly. Without them the inlet takes
I = 5 % and ν_t/ν = 10 (TURB009 names it). Rough walls take the
sand-grain roughness:

```toml
[boundary.inlet]
type = "velocity-inlet"
velocity = [10.0, 0.0, 0.0]
turbulence_intensity = 0.02
turbulence_length_scale = 0.01   # or: viscosity_ratio = 10.0, or: k = ..., omega = ...
[boundary.wall]
type = "wall"
roughness = 1e-4                 # k_s (m), optional
```

The start values come from the first inlet unless `[initial]` gives
`k` and `omega` (numbers or expressions) or `turbulence_intensity` and
`viscosity_ratio`. Check the case, run it and read the first-cell y⁺ per
wall, the largest ν_t/ν and the floor counts at the end of the log:

```
target/release/foehn flow examples/cases/flow/tmr_plate_sst.toml --check
target/release/foehn flow examples/cases/flow/tmr_plate_sst.toml
grep -A3 "^turbulence:" output/tmr_plate_sst/flow.log
```

The report has the same under `turbulence` in `result.json`; the wall
shear stress and y⁺ per wall face are in `wall_<boundary>.csv`
(`output.wall_csv = true`). For a mesh study of C_f, run the case on four
meshes and read the discretisation uncertainty in `flow.md`:

```
target/release/foehn flow examples/cases/flow/tmr_plate_sst.toml --out output/plate \
    --set 'verification.levels=[0, 1, 2, 3]'
```

A fully developed channel needs one periodic module (recipe 67) and the
`low-re` treatment with y⁺ ≲ 1:

```
target/release/foehn flow examples/cases/flow/channel_sst.toml
grep "periodic xmin" output/channel_sst/flow.log
```

u_τ = √(β h/ρ) from the driving force β and the half height h, so
Re_τ = u_τ h/ν and U_b⁺ = U_b/u_τ. `output.derived = ["turbulence-features"]`
writes the invariant features of D166 (the inputs of data-driven
corrections) with the fields. The wall distance is exact by default;
`turbulence.wall_distance = "pde"` takes the Poisson distance instead
(cheaper on very large meshes, exact on plane walls).

[↑ back to top](#turbulence-models)

---

<a id="k-epsilon"></a>
## 71. Run a case with the k-ε family and wall functions

The standard k-ε (`k-epsilon`) or the realizable k-ε
(`realizable-k-epsilon`) with wall functions: put the first cell in the
log layer, y⁺ ≳ 30. `scalable` (the default) takes any y⁺ without
failing and stays close below 30; `standard` needs 30 ≤ y⁺ ≤ 300;
`compound` is offered, but with these high-Re models it is claimed from
y⁺ 30 only. Walls resolved to y⁺ ≈ 1 need SST `low-re` now, ζ-f or EB
later (TURB001 warns):

```toml
[turbulence]
model = "k-epsilon"              # or "realizable-k-epsilon"
wall_treatment = "scalable"      # "standard", "compound"
# production_limiter = "durbin"  # the default of the standard model; "kato-launder", "menter" (min(P, 10ε)), "none"
# ambient = { k = 2.25e-7, epsilon = 2.53125e-6 }   # sustaining terms, external flows
```

Inlets take ε directly or ω (converted with ε = C_μ k ω), or the
intensity with a length scale or viscosity ratio, as for SST:

```toml
[boundary.inlet]
type = "velocity-inlet"
velocity = [1.0, 0.0, 0.0]
k = 2.25e-7
epsilon = 2.53125e-6             # or: omega = 125.0
```

With free-stream values that decay fast (the TMR plate: ν_t/ν = 0.009)
give the sustaining terms the inlet values, or the boundary layer starts
from nothing. Run the plate, the channel and the rough pipe of the
gallery:

```
target/release/foehn flow examples/cases/flow/tmr_plate_keps.toml
target/release/foehn flow examples/cases/flow/channel_keps.toml
target/release/foehn flow examples/cases/flow/pipe_rough_keps.toml
grep -E "wall probe|periodic xmin|first-cell" output/*_keps/flow.log
```

C_f of the plate is in the wall probe line; in the channel
U_b⁺ = U_b/√(F h/ρ) with the driving force F; in the pipe the Darcy
friction factor is f = 2 F D/(ρ U²). A run that comes from SST starts
from its restart file: `[run] restart = 'output/plate/restart.h5'` with
`model = "k-epsilon"` converts ω to ε.

Rough walls: `roughness` is the sand-grain k_s; C_s of Cebeci & Bradshaw
is `roughness_constant` (0.5). For pipes with Colebrook-type commercial
roughness, 0.4 puts the fully rough friction factor within 3 % of
Colebrook's (0.5: 8 to 11 % above):

```toml
[boundary.wall]
type = "wall"
roughness = 1e-3
roughness_constant = 0.4
```

The Moody chart as a parametric study (16 runs, about 8 minutes with four
jobs; needs the Python package or `FOEHN=target/release/foehn`):

```
python3 examples/studies/moody_study.py output/moody --jobs 4 --chart output/moody.svg
python3 examples/studies/moody_study.py output/moody_cs04 --cs 0.4 --jobs 4
```

The y⁺ study of the plate and the figure of gallery 71:

```
for r in "1 17193.2" "2 7983.56" "5 2859.14" "10 1298.87" "30 361.051" "60 156.835" "100 83.3061"; do set -- $r
  for t in scalable compound; do
    target/release/foehn flow examples/cases/flow/tmr_plate_keps.toml --out output/plate_keps/${t}_$1 \
        --set "turbulence.wall_treatment=\"$t\"" --set "mesh.generate.grading.y.ratio=$2"
  done
done
python3 -m tools.plot.keps_chart docs/figures/gallery/keps.svg output/plate_keps output/moody/study.csv
```

[↑ back to top](#turbulence-models)

---

<a id="zeta-f"></a>
## 73. Run a case with k-ε-ζ-f and resolved walls

k-ε-ζ-f (`k-epsilon-zeta-f`) for walls resolved to y⁺ ≈ 1, where the
k-ε family needs wall functions. `compound` (the default) takes any y⁺;
`low-re` needs y⁺ ≲ 2 (TURB001):

```toml
[turbulence]
model = "k-epsilon-zeta-f"
wall_treatment = "compound"      # or "low-re"
# ambient = { k = 2.25e-7, epsilon = 2.53125e-6 }   # sustaining terms, external flows
```

Inlets take k with ε (or ω, converted) as for k-ε, and ζ = v²/k
(default 2/3, isotropic turbulence); `[initial]` takes `zeta` too:

```toml
[boundary.inlet]
type = "velocity-inlet"
velocity = [1.0, 0.0, 0.0]
k = 1.0e-4
epsilon = 4.5e-5
zeta = 0.667
```

A weak free stream (ν_t/ν well below 1) leaves the boundary layer
laminar with this model: start from a turbulent field (`[initial] k`,
`epsilon` with ν_t/ν ≈ 10), as the plate of gallery 73 does, or give the
inlet more turbulence.

The f equation runs CG with its own multigrid hierarchy; the transport
equations take the multigrid too, which pays on large meshes:

```toml
[solver.elliptic]
preconditioner = "amg"           # default; "jacobi"
tolerance = 1e-2                 # default steady (transient 1e-8)

[solver.turbulence]
preconditioner = "amg"           # default "jacobi"; with the coupled algorithm
                                 # "reduction" (recipe 81)
```

Run the channel and the plate of the gallery:

```
target/release/foehn flow examples/cases/flow/channel_zeta_f.toml
target/release/foehn flow examples/cases/flow/tmr_plate_zeta_f.toml
grep -E "wall probe|periodic xmin|first-cell|floors" output/*_zeta_f/flow.log
```

U_b⁺ of the channel is 17.55/√F with the driving force F; C_f of the plate
is in the wall probe line. The profile line of the channel
(`line_profile.csv`) has ζ and f next to k, ε and ν_t; the fields add
`zeta`, `f`, `turbulent_time_scale` and v² = ζk. A k-ε restart file
continues with ζ = 2/3 and f = 0.

[↑ back to top](#turbulence-models)

---

<a id="zeta-f0"></a>
## 74. Run a case with k-ε-ζ-f₀ at any first-cell y⁺

k-ε-ζ-f₀ (`k-epsilon-zeta-f0`, Popovac 2023) for meshes whose first-cell
y⁺ is not controlled, from resolved walls to the log layer.
`compound-2023` is its default treatment; `compound` and `low-re` (y⁺ ≲ 2)
are also accepted, and ζ-f takes `compound-2023` too:

```toml
[turbulence]
model = "k-epsilon-zeta-f0"
wall_treatment = "compound-2023"   # default; "compound", "low-re"
# k_sink = "wall"   # the default (recipe 141)
# constants = { n_zeta = 2.0, c_epsilon1_zeta = 0.009, c_epsilon1_power = 1.0 }   # the defaults
```

Inlets and `[initial]` take k, ε (or ω) and ζ as for ζ-f (recipe 73).
Unlike ζ-f, a weak free stream is fine: the TMR plate starts from
ν_t/ν = 0.009 and turns turbulent on its own.

A step or another geometry with a cut-out corner can come from the
built-in mesher: `remove` cuts boxes out of the mesh and names the new
faces:

```toml
[mesh.generate]
type = "hex"
dimension = 2
lower = [-85.0, 0.0, 0.0]
upper = [40.0, 9.0, 0.1]
cells = [250, 160, 1]
remove = [ { lower = [-85.0, 0.0, 0.0], upper = [0.0, 1.0, 0.1], name = "step" } ]

[boundary.step]
type = "wall"
```

Run the channel, the plate and the step of the gallery:

```
target/release/foehn flow examples/cases/flow/channel_zeta_f0.toml
target/release/foehn flow examples/cases/flow/tmr_plate_zeta_f0.toml
target/release/foehn flow examples/cases/flow/backstep_zeta_f0.toml
grep -E "wall probe|periodic xmin|first-cell|floors" output/*_zeta_f0/flow.log
```

The reattachment length of the step is where the wall shear of `ymin`
in `wall_ymin.csv` (behind the step, x > 0) changes sign from negative
to positive.

Compare the models on the step. The three cases in
`examples/studies/backstep/` are the step on the fine mesh (32 800
cells) with k-ε and scalable wall functions, k-ε-ζ-f₀ and k-ω SST, each
with `output.vtkhdf = true`:

```
target/release/foehn flow examples/studies/backstep/keps.toml --out output/backstep_keps
target/release/foehn flow examples/studies/backstep/zf0.toml --out output/backstep_zeta_f0
target/release/foehn flow examples/studies/backstep/sst.toml --out output/backstep_sst
```

Each case starts with the momentum under-relaxation 0.5 while the bubble
forms and switches to 0.8 at iteration 100:

```toml
[solver.momentum]
relaxation = 0.5

[[schedule]]
at_iteration = 100
set = { solver.momentum.relaxation = 0.8 }
```

That halves the iterations: x_r is within 1 % of its final value after
about 600 iterations on the 8 200-cell mesh and 1 750 on the 32 800-cell
mesh, for all three models. With 0.5 throughout SST needs more than
4 000 on the fine mesh; 0.9 stalls the residual at 10⁻⁴. The residual
is not the measure here: it is still 5·10⁻⁵ to 10⁻⁴ when x_r has
settled, so the cases run to `max_iterations` (3 000) with
`require_convergence = false`. To see when a result settles, have the
run write the fields at a few iterations and compare them:

```
mkdir -p output/backstep_sst
for n in 250 500 1000 1500 2000; do touch output/backstep_sst/.result_$n; done
target/release/foehn flow examples/studies/backstep/sst.toml --out output/backstep_sst
```

The figures of the gallery, velocity with streamlines and turbulent
kinetic energy, one panel per model, come from ParaView's Python or from
matplotlib. Both scripts take an output prefix and `LABEL=DIR` per
model, write `PREFIX_velocity.png` and `PREFIX_k.png` and print x_r/H:

```
pvpython tools/plot/backstep_paraview.py docs/figures/gallery/backstep \
    "k-ε (scalable)=output/backstep_keps" "k-ε-ζ-f₀=output/backstep_zeta_f0" "k-ω SST=output/backstep_sst"
python3 -m tools.plot.backstep_chart docs/figures/gallery/backstep \
    "k-ε (scalable)=output/backstep_keps" "k-ε-ζ-f₀=output/backstep_zeta_f0" "k-ω SST=output/backstep_sst"
```

On macOS `pvpython` is in the application,
`/Applications/ParaView-6.2.0.app/Contents/bin/pvpython`. The ParaView
script is a starting point for other 2D cases: it reads `fields.vtkhdf`,
cuts the mesh at mid-depth, and sets the camera window, the colour map
and the axis text itself, so the picture does not depend on the ParaView
version (tested with 5.11 and 6.2).

[↑ back to top](#turbulence-models)

<a id="elliptic-blending"></a>
## 88. Run a case with the elliptic-blending k-ε

The elliptic-blending k-ε of Billard & Laurence (2012)
(`elliptic-blending-k-epsilon`) resolves the wall like ζ-f: keep the
first-cell y⁺ near 1 or below. `compound` is its default treatment,
`low-re` (y⁺ ≲ 1) the alternative; the wall functions of the k-ε family
and `compound-2023` are refused (CASE005):

```toml
[turbulence]
model = "elliptic-blending-k-epsilon"
wall_treatment = "compound"   # default; "low-re"
```

Inlets and `[initial]` take k and ε (or ω) as for ζ-f; the key `zeta`
sets φ = v²/k (default 2/3). α needs no value: it starts at 1 and
follows from k and ε. The log shows the residuals `ResPhi` and
`ResAlpha`; the output fields are `phi` and `alpha`.

Low Reynolds numbers need a strong start: the channel at Re_τ 180 turns
laminar from I = 0.05 and ν_t/ν = 100 and stays turbulent from

```toml
[initial]
turbulence_intensity = 0.15
viscosity_ratio = 20.0
```

Separated flows: the paper's variable C*_ε2 lengthens the bubble of the
backward-facing step to x_r/H = 7.59 (+21 %). A constant C_ε2 gives
6.24, at the price of +4.9 % on the C_f of the flat plate:

```toml
[turbulence]
constants = { c_epsilon4 = 1.83 }   # C*_ε2 = C_ε2 = 1.83 everywhere
```

On the step, use SIMPLEC with the relaxation ramp of
`backstep_eb.toml`; with `algorithm = "coupled"` the step does not
converge (the residual of ε stays at 10⁻³).

Run the channel, the plate and the step of the gallery:

```
target/release/foehn flow examples/cases/flow/channel_eb.toml
target/release/foehn flow examples/cases/flow/tmr_plate_eb.toml
target/release/foehn flow examples/cases/flow/backstep_eb.toml
grep -E "wall probe|periodic xmin|first-cell|floors" output/*_eb/flow.log
```

The channels at Re_τ 590 (128 cells, grading 80, first-cell y⁺ 0.25)
and 180 of the gallery, with ν = 1/Re_τ and the mass flow of the DNS
bulk velocity; C_f/C_f,DNS is the driving force in the log:

```
target/release/foehn flow examples/cases/flow/channel_eb.toml --out output/channel_eb_590 \
    --set 'mesh.generate.cells=[2,128,1]' --set mesh.generate.grading.y.ratio=80.0 \
    --set fluid.viscosity=0.001694915254237288 --set boundary.xmin.mass_flow=3.7124
target/release/foehn flow examples/cases/flow/channel_eb.toml --out output/channel_eb_180 \
    --set fluid.viscosity=0.005555555555555556 --set boundary.xmin.mass_flow=3.102 \
    --set initial.turbulence_intensity=0.15 --set initial.viscosity_ratio=20.0
grep "periodic xmin" output/channel_eb_590/flow.log output/channel_eb_180/flow.log
```

[↑ back to top](#turbulence-models)

<a id="step-validation"></a>
## 89. Validate the turbulence models on the backward-facing step

The validation case of stage 6 (specification 9.5) is the step of
Driver & Seegmiller (1985) in the layout of NASA's Turbulence Modeling
Resource: `backstep_validation.toml`, k-ε-ζ-f₀ with the coupled
algorithm, no-slip walls at the bottom and the top, 11 040 cells. One
run of it:

```
target/release/foehn flow examples/cases/flow/backstep_validation.toml
grep -E "wall crossing|wall probe" output/backstep_validation/flow.log
```

**The reattachment length as a result.** A wall crossing reports the
first point after `from` along `direction` where the wall shear along
`direction` changes sign, linear between the two face centroids; its
distance from `from` is in the log, in `flow.json` and, in a sweep over
three or more meshes, in the uncertainty table with the observed order
and the GCI (D67):

```toml
[output]
wall_crossings = [
  # reattachment (the default kind), searched from x = 2 H on
  { name = "xr", boundary = "ymin", from = [0.0, 0.0, 0.05], skip = 2.0 },
  { name = "xs", boundary = "ymin", from = [0.0, 0.0, 0.05], kind = "separation" },
]
```

`direction` defaults to x; `skip` leaves out the faces closer than that
distance to `from` and still measures the position from `from`: behind
the step the corner vortices turn the shear positive and back up to
x ≈ 1.5 H on the finer meshes, and without `skip` the first crossing is
theirs. `width` keeps only the faces within that distance of the line
through `from` (for a 3D wall). The faces are taken
in the order of their centroids along the direction, so the wall should
be one row of faces: a 2D case, or a 3D case with `width`.

**Three meshes.** `verification.levels` halves every interval once and
twice (44 160 and 176 640 cells); the uncertainty table of `flow.md`
then has the row `wall_crossings.xr.position`:

```
target/release/foehn flow examples/cases/flow/backstep_validation.toml --out output/step_gci \
    --set 'verification.levels=[0, 1, 2]' --set solver.max_iterations=800
grep -A12 "Discretisation uncertainty" output/step_gci/flow.md
```

**All models against the experiment.** The study runs SST, k-ε,
realizable k-ε, ζ-f, ζ-f₀ and the elliptic-blending k-ε on the three
meshes (ζ-f and EB with SIMPLEC and the relaxation ramp; the coupled
algorithm does not converge them), then compares the finest mesh with
the experiment: x_r/H with its GCI, the RMS difference of C_p at the 33
stations of the lower wall (the band of the specification is 0.02 for
SST, ζ-f and EB), C_f, and the velocity profiles at x/H = 1, 4, 6 and
10. The reference speed is u at x/H = −4, y/H = 5, and C_p is shifted to
zero at x/H = 37.5 like the data:

```
pip install ./python
python3 examples/studies/backstep_validation.py output/backstep_validation --jobs 4 \
    --chart docs/figures/gallery
cat output/backstep_validation/summary.md
```

`--models sst,keps,rke,zf0` runs a subset; `--levels 0,1` fewer meshes;
`--report-only` evaluates the runs that are there. A run whose case has
not changed is not repeated, so the study can be stopped and started
again. On four cores the four coupled models take about an hour
together (k-ε 23 minutes for its three meshes, SST 36); ζ-f and EB need
hours with SIMPLEC on the finest mesh, and `--models zf,eb --levels 0,1`
runs their two coarser meshes in about an hour.

[↑ back to top](#turbulence-models)

<a id="impinging-jet"></a>
## 90. Feed a jet with a developed pipe profile and compute the Nusselt number

The jet of specification 9.6 leaves a long pipe, so its inlet carries
the developed velocity *and* turbulence of a pipe flow. Two cases do it:
a periodic pipe on a wedge (`jet_pipe_profile.toml`, 0.3 s) and the jet
on a heated plate (`impinging_jet.toml`, about a minute on four cores).

**1. The developed profile.** The periodic pipe at the jet's Reynolds
number samples the radius with a line from the axis to the wall:

```
target/release/foehn flow examples/cases/flow/jet_pipe_profile.toml --out output/jet_pipe
```

**2. The inlet table.** The line samples the turbulence variables as
cell values; `foehncraft.profiles` takes one row per cell (the centre of
its run of equal samples, the velocity there) and adds a wall row
(u = 0, k = 0, ζ = 0, ε of the wall cell):

```
python3 -m foehncraft.profiles output/jet_pipe/line_profile.csv \
    examples/cases/flow/impinging_jet_profile.csv
```

The header of the table names its columns: the coordinate, then the
velocity as `u_n` (along the inlet normal; here along
`table_direction`) or `u_x`, `u_y`, `u_z`, and the turbulence as `k`
with `epsilon` or `omega`, and `zeta`:

```toml
[boundary.xmax]
type = "velocity-inlet"
profile = "table"
table = "impinging_jet_profile.csv"   # r, u_n, k, epsilon, zeta
table_coordinate = "radius"
table_direction = [1.0, 0.0, 0.0]
```

With turbulence columns the boundary gives no turbulence values of its
own; `turbulence_intensity`, `k`, … or `zeta` next to such a table is
refused (BC029), and so is a column name the profile does not know
(BC019). A table of another model works the same way (`omega` for SST;
`zeta` is φ for the elliptic-blending k-ε).

**3. The Nusselt number.** The plate is a wall with a uniform flux of
the passive temperature θ (Pr = 0.71), and both inlets bring θ = 0, so
Nu(r) = q D/(Γ θ_w(r)) with Γ = ν/Pr from the column `theta_wall` of
`wall_xmin.csv`. A wall probe reports the wall value and the flux of
every scalar at a point, here at the stagnation point and at r/D = 0.1:

```toml
[output]
wall_csv = true
wall_probes = [
  { name = "stag", boundary = "xmin", point = [0.0, 0.0, 0.0] },
  { name = "r01", boundary = "xmin", point = [0.0, 0.1, 0.0] },
]
```

```
target/release/foehn flow examples/cases/flow/impinging_jet.toml
grep "wall probe" output/impinging_jet/flow.log
```

The log prints both (`wall probe stag at (0, 2.94e-2, 0) m: … theta
wall 1.075e-1, flux 1.000e-3`). A probe interpolates between the two
wall faces nearest its point and is clamped to the nearest face beyond
them, so `stag` reports the face next to the axis (r/D = 0.03 on this
mesh). In a sweep over three meshes `wall_probes.stag.theta_wall` and
`wall_probes.r01.theta_wall` get their observed order and GCI (D67), and
so Nu₀ (to first order the same relative GCI). The two differ by 0.5 %
on this mesh and by less than 0.1 % on 41 600 cells. (Before PR #91 the
faces next to the axis were up to 7 % low with the turbulence models:
the planes of the wedge lacked the transposed part of the stress, half
of the hoop stress; see the verification page.)

**All models, three meshes, the limiters.** The study runs the pipe and
the jet for six models (`verification.levels` 0, 1, 2: 2 600, 10 400
and 41 600 cells) and the cases of the stagnation-point anomaly (k-ε and
SST with each production limiter on the 10 400 cells), then writes
`summary.md`; the experiment (ERCOFTAC case 025) enters when its CSV
files are in `examples/studies/impinging_jet/data` (see the README there;
`ercoftac_to_csv.py` makes them), and the charts then draw it, so the
gallery figures are made with `--no-experiment`:

```
pip install ./python
python3 examples/studies/impinging_jet.py output/impinging_jet_study --jobs 2 \
    --chart docs/figures/gallery
cat output/impinging_jet_study/summary.md
```

`--models zf0,sst` runs a subset, `--levels 1,2,3` the finer meshes
(166 400 cells on the finest), `--no-limiters` leaves the anomaly out
and `--report-only` evaluates what is there. ζ-f and EB run with SIMPLEC
and the relaxation ramp of the step (the coupled algorithm does not
converge them here either).

[↑ back to top](#turbulence-models)

---

<a id="curvature-correction"></a>
## 113. Make SST sensitive to curvature and system rotation

Linear eddy-viscosity models do not feel streamline curvature or system
rotation: in a swirling, rotating or strongly curved flow they keep the
turbulence of a plain shear layer. The correction of Smirnov & Menter
(2009) multiplies the production of k and ω by f_r1 (0 … 1.25):

```toml
[turbulence]
model = "k-omega-sst"
curvature_correction = true
# constants = { c_r1 = 1.0, c_r2 = 2.0, c_r3 = 1.0 }   # the paper's values
```

- It acts in the fixed frame (curvature) and sees the rotation of
  `[[mrf]]` zones (the velocity relative to the zone convects the strain
  rate, Ω of the zone enters r̃). It does not work with `[[sliding]]` yet.
- A flow in system rotation that is periodic in the rotating frame (a
  rotating channel or duct) is solved in that frame with
  `flow.frame_rotation = [Ωx, Ωy, Ωz]` (rad/s): the velocity is the
  relative one and the momentum equations get the Coriolis force; the
  centrifugal force is a gradient and stays in the pressure.
- `output.derived = ["curvature-correction"]` (or `foehn post derive
  --fields curvature-correction`) writes f_r1 per cell: blue (< 1) where
  the flow damps the turbulence, red (> 1) where it amplifies it.

```
target/release/foehn flow examples/cases/flow/channel_rotating_sst.toml
```

The run logs the rotating frame; `flow.md` lists the wall forces of the
two walls (`ymin` is the unstable side at Ω_z > 0). Without the
correction both walls carry the same shear stress; with it they differ by
a factor 1.63 at Ro = 0.15. The cost is six extra gradients per outer
iteration. Combined with `production_limiter = "kato-launder"` the two
act on the same regions (stagnation and swirl); the paper uses the
standard production.

[↑ back to top](#turbulence-models)

---

<a id="qcr"></a>
## 115. Get the secondary flow of non-circular ducts with SST

A linear eddy-viscosity model gives no secondary flow in a straight
square or rectangular duct, nor in the corners of blades and walls. The
quadratic constitutive relation makes the normal stresses anisotropic,
which drives the corner vortices:

```toml
[turbulence]
model = "k-omega-sst"
qcr = "qcr2000"          # or "qcr2013"; "none" is the default
# constants = { c_cr1 = 0.3, c_cr2 = 2.5 }   # the papers' values
```

- The non-linear stress is explicit. It is added to the momentum
  equations from the velocity gradient of the current iteration and ν_t
  of the last turbulence update. The implicit diffusion stays μ + μ_t, and
  the production of k keeps the linear form, as in the papers.
- `"qcr2013"` adds an isotropic stress. In incompressible flow that only
  changes the pressure, so the velocity is the one of QCR2000.
- For SST only. It combines with `curvature_correction`.
- The secondary velocity is about 1.5 % of U_b, so converge the
  in-plane residuals (V, W) well. In the gallery duct they fall like the
  others and the run needs about 2 600 iterations on 10 368 cells.

```
target/release/foehn flow examples/cases/flow/square_duct_qcr.toml
```

[↑ back to top](#turbulence-models)

---

<a id="k-sink"></a>
## 141. Choose the sink of k of the ζ-f models

k-ε-ζ-f and k-ε-ζ-f₀ take the sink max(k/T_K, 2νk/d²) of k by default
(`"wall"`, d the wall distance): it is ε at the wall, so the results
converge as the first cell moves towards the wall, and k/T_K away from
it, so k stays positive in a laminar boundary layer. Nothing needs to be
set. The other forms are there to compare with the literature and with
older results:

```toml
[turbulence]
model = "k-epsilon-zeta-f0"
k_sink = "epsilon"       # −ε as printed (Hanjalić et al. 2004, Popovac 2023)
# k_sink = "kolmogorov"  # k/T_K, the ζ-f₀ default before #141
```

The papers' −ε holds k at its floor at a laminar leading edge (the k
residual stalls, TURB002 counts the floors in the log); `"kolmogorov"`
lets k fall only like y at the wall, so U_b⁺ of the channel drifts with
the mesh below a first-cell y⁺ of 1. ζ-f₀ before #141 is

```toml
[turbulence]
model = "k-epsilon-zeta-f0"
k_sink = "kolmogorov"
constants = { c_epsilon1_zeta = 0.012, c_epsilon1_power = 0.5 }   # 1.4 (1 + 0.012/√ζ)
```

Check that a result converges with the mesh with the study of gallery
141 (the channel refined in y on four levels, the coupled algorithm):

```
python3 examples/studies/plate_channel_gci.py output/gci --cases channel --models zf0,zf --levels 0,1,2,3
```

EB (`elliptic-blending-k-epsilon`) keeps −ε and refuses `k_sink`
(CASE005).
<a id="outlet-backflow"></a>
## 140. Say what turbulence enters through an outlet

When fluid enters through part of a pressure outlet (an outlet close
behind a body, in a recirculation, at an open side of a domain), give the
turbulence of the entering fluid with the keys of an inlet:

```toml
[boundary.outlet]
type = "static-pressure-outlet"      # or "total-pressure-outlet"
pressure = 0.0
turbulence_intensity = 0.01          # or k with omega / epsilon
viscosity_ratio = 10.0               # or turbulence_length_scale
# zeta = 0.6667                      # ζ-f, ζ-f₀, EB (default 2/3)
```

- Only the faces where the flow enters (ṁ < 0) take them, as fixed values;
  the others stay zero gradient. The split is made anew in every outer
  iteration.
- `turbulence_intensity` takes the speed of the face, like at an inlet.
- Without values the entering fluid carries the turbulence of the cell
  next to the outlet, which copies a shear layer into the inflow.
- The momentum treatment of the same faces is the separate key
  `backflow` ([recipe 29](boundary-conditions.md#boundary-conditions)).
  `backflow = "suppress"` closes the faces, and then the turbulence values
  are not used.
- Better still, move the outlet further downstream; the values only keep
  an unavoidable inflow plausible.

```
target/release/foehn flow examples/cases/flow/backstep_outlet_backflow.toml
```

[↑ back to top](#turbulence-models)

---

<a id="qcr-consistency"></a>
## 194. Run the QCR on tetrahedra and perturbed meshes

Nothing to switch on: a run with `turbulence.qcr` checks its velocity
gradient before the solve (D236) and uses the quadratic least-squares
gradient for the QCR stress on every cell where the face-neighbour
gradient is first order. The log says how many cells that is, and
`foehn flow --check --json` shows it before a run:

```
target/release/foehn flow examples/cases/flow/cavity_sst_operator_check.toml --check --json
```

- `operators.quadratic` in the report holds the error of the quadratic
  fit; `MESH006` remains only for the cells where the fit fails (too few
  or degenerate neighbours of neighbours). Remesh those cells, or convert
  tetrahedra to polyhedra (D54).
- `solver.operator_tolerance` (0.01) sets which cells count as flagged and
  so which cells take the quadratic gradient.
- On regular hexahedra the quadratic gradient is used only in the cell
  layer at the walls; inside the mesh the gradient is already exact.
- The curvature correction (`turbulence.curvature_correction`) takes S_ij
  from the same fit (recipe 197).

```
target/release/foehn flow examples/cases/flow/cavity_sst_operator_check.toml
```
<a id="pipe-turbulent"></a>
## 193. Run a turbulent pipe and check its friction

Use the wedge pipe (recipe 71) with a periodic module and `low-re`:

```
target/release/foehn flow examples/cases/flow/pipe_turbulent_sst.toml --out output/pipe
grep "periodic xmin" output/pipe/flow.log
```

With the body force F of the log, u_τ = √(F R/2ρ) (ρ = 1, R = 0.5 m),
u_τ/U_b with U_b = 1 and C_f = 2 (u_τ/U_b)²; Re_τ = u_τ R/ν. The case
gives u_τ/U_b of 0.0711 on 96 cells (DNS 0.0683). For the grid-convergence
study (observed order, extrapolated value, GCI_fine after Celik et al.
2008) run `python3 examples/studies/pipe_gci.py output/pipe_gci`; it
refines across the radius only, the solution does not depend on x. Steady
wall-resolved runs need `solver.momentum.relaxation = 0.95` here, and
judge convergence by u_τ, not by the residual alone.

[↑ back to top](#turbulence-models)

---

<a id="diffuser"></a>
## 199. Run an asymmetric diffuser, its grid study and a transient run from the steady state

The inlet is a developed channel profile: the periodic channel
`diffuser_channel_profile.toml` (bulk velocity set so that the centreline
velocity is 1 m/s, Re_H = U_c H/ν = 20 000) gives the table
`diffuser_inlet_profile.csv` (`python3 examples/studies/diffuser_inlet_profile.py`;
change the mass flow of the periodic pair if the viscosity changes).
Make the structured mesh at a level (every level halves all intervals;
both corners of the lower wall are rounded with 4.3 H) and run:

```
python3 examples/meshes/diffuser.py diffuser.su2 --level 0
target/release/foehn flow examples/cases/flow/diffuser_sst.toml --out output/diffuser \
    --set 'mesh.file="diffuser.su2"'
```

- The case starts from potential flow (`[initial] potential_flow`), not
  from a uniform field, with first-order upwind momentum, and switches to
  linear upwind at iteration 201:

  ```toml
  [solver.momentum]
  scheme = "upwind"
  relaxation = 0.8

  [[schedule]]
  at_iteration = 201
  set = { solver.momentum.scheme = "linear-upwind" }
  ```

  The scheme may change only between `"upwind"` and `"linear-upwind"`
  (also at run time through `control.toml`); a change to or from
  `"central"` is refused.
- The report prints the separation and the reattachment point of the
  lower wall (`output.wall_crossings`, positions from x = 0) and the
  first-cell y⁺ of both walls; `wall_upper.csv` holds the wall pressure
  for the pressure recovery. The log prints u at a probe in the jet core
  and one in the separated region: judge convergence on these. On the
  finer levels the momentum residual rises for a while (on level 2 from 1 000 to 2 400 of its 4 909 iterations) while
  the jet moves to the upper wall; it is not a divergence.
- **Sharp corners**: `diffuser.py --sharp` makes the mesh without the
  radius and shortens the diffuser cells towards both corners (aspect
  ratio about 11 at the kinks); run it with `diffuser_sst_sharp.toml`. Its
  separation search starts 1 H behind the corner (`skip = 1.0`), since the
  wall shear stays near zero there and small bubbles come and go with the
  mesh. For a sharp corner of your own geometry, keep the first-order
  start and shorten the wall cells at the kink
  ([pitfalls](../pitfalls.md#sharp-wall-corners-and-thin-wall-cells-the-asymmetric-diffuser)).
- The grid study runs levels 0 to 2 at the same time and writes
  `report.md` with the order and GCI of each quantity:

  ```
  python3 examples/studies/diffuser_gci.py output/diffuser_gci --jobs 3
  python3 examples/studies/diffuser_gci.py output/diffuser_gci_sharp --jobs 3 --sharp
  ```

- A **transient run from the steady state**: give the steady restart file
  and a time step. Drop the iteration schedule (`at_iteration` does not fit
  a transient case) and set the relaxation to 1 (one outer iteration per
  step):

  ```
  target/release/foehn flow examples/cases/flow/diffuser_sst.toml --out output/diffuser_t \
      --set 'mesh.file="diffuser.su2"' --restart output/diffuser/restart.h5 \
      --set 'schedule=[]' --set solver.momentum.relaxation=1.0 \
      --set time.time_step=0.1 --set time.end_time=100.0 \
      --set 'output.probes=[{name="sep",point=[10.0,-1.5,0.5]},{name="shear",point=[15.0,-1.0,0.5]},{name="reat",point=[28.0,-3.5,0.5]},{name="core",point=[15.0,0.5,0.5]}]'
  python3 -m tools.plot.diffuser_transient transient.png output/diffuser_t OTHER_RUN
  ```

  The log notes "start from the steady solution of …"; the steady fields
  are the start fields at `time.start_time`, with no older time levels.
  `probes.csv` holds the four probes; a steady state that is stable stays
  flat (gallery 199).

<a id="curvature-consistency"></a>
## 197. Run the curvature correction on tetrahedra and perturbed meshes

Nothing to switch on: with `turbulence.curvature_correction = true` the
run checks the gradient of its velocity gradient before the solve (D236)
and, on every cell where it is not exact for a quadratic field, takes the
gradient of S_ij in DS_ij/Dt from a quadratic least-squares fit of S_ij
(with S_ij from the quadratic velocity fit of recipe 194 where the
face-neighbour gradient is first order). The log names the number of
cells:

```
target/release/foehn flow examples/cases/flow/cavity_sst_operator_check.toml
# curvature correction: S_ij from the quadratic velocity fit on 4096 cells, DS_ij/Dt from a quadratic least-squares fit of S_ij on 4096 cells (D236)
```

- `foehn flow --check --json` shows it before a run: `operators.hessian`
  is the check of the plain operator, `operators.quadratic_hessian` that of
  fits on the treated cells. `MESH005` remains only for cells where a
  fit fails (too few or degenerate neighbours of neighbours); remesh those
  cells or convert tetrahedra to polyhedra (D54).
- On regular hexahedra and polyhedra only the cell layers at the walls
  take the fit; the result changes only there.
- With the QCR on as well, both use the same velocity fit.

[↑ back to top](#turbulence-models)

<a id="ubend"></a>
## 204. Run a U-bend, its grid study and the sharp-edged variant

The inlet is a developed channel profile: the periodic channel
`ubend_channel_profile.toml` (mass flow of U_b = 1 m/s, Re_b = 10⁵) gives
the table `ubend_inlet_profile.csv`, already placed on the inlet leg at
−1.5 < y < −0.5 (`python3 examples/studies/ubend_inlet_profile.py`; run it
again if the viscosity changes). Make the structured mesh at a level (every
level halves all intervals) and run:

```
python3 examples/meshes/ubend.py ubend.su2 --level 0
target/release/foehn flow examples/cases/flow/ubend_sst.toml --out output/ubend \
    --set 'mesh.file="ubend.su2"'
```

- The case starts from potential flow with first-order upwind momentum and
  switches to linear upwind at iteration 201, relaxation 0.8 (as the
  diffuser, recipe 199).
- The report prints the reattachment point on the inner wall behind the
  turn (`output.wall_crossings` from the end of the turn in −x; 2.40 H on
  level 0, within one wall face of the 2.32 H the study interpolates from
  the wall table) and the first-cell y⁺; `wall_inner.csv` and `wall_outer.csv` hold the wall
  pressure and shear. The log prints u at a probe in the core of the outlet
  leg and one in the separated region (u > 0 there is reverse flow, the
  outlet leg runs in −x): judge convergence on these.
- The grid study runs levels 0 to 2 at the same time and writes
  `report.md` with the order and GCI of the separation angle, the
  reattachment and three wall pressures. It measures positions along the
  wall (the angle in the turn, the distance behind it) and takes the
  shear along the flow direction, so a separation inside the turn is found
  too:

  ```
  python3 examples/studies/ubend_gci.py output/ubend_gci --jobs 3
  python3 examples/studies/ubend_gci.py output/ubend_gci_sharp --jobs 3 --sharp
  ```

- **Sharp edges**: `ubend.py --sharp` makes a square divider end and a
  square end wall with square cells at the kinks; run it with
  `ubend_sst_sharp.toml`. It starts at relaxation 0.5, switches to linear
  upwind at iteration 601 and raises the relaxation to 0.8 at 1 201:

  ```toml
  [solver.momentum]
  scheme = "upwind"
  relaxation = 0.5

  [[schedule]]
  at_iteration = 601
  set = { solver.momentum.scheme = "linear-upwind" }

  [[schedule]]
  at_iteration = 1201
  set = { solver.momentum.relaxation = 0.8 }
  ```

  0.8 from the start diverges at the concave corners, and switching both at
  once diverges behind the divider end; 0.5 throughout converges too
  slowly on level 1.
- The curvature correction of SST changes the separation behind the turn
  (gallery 204): `--set turbulence.curvature_correction=true`.
- Figures: `python3 -m tools.plot.ubend_chart OUT.png [WORKDIR] [--sharp]`
  (fields and wall pressure of the coarse mesh) and
  `python3 -m tools.plot.ubend_friction OUT.png GCI_DIR CC_DIR SHARP_DIR`
  (inner-wall friction of the study runs).

[↑ back to top](#turbulence-models)
