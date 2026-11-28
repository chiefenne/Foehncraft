# Lagrangian sprays

Part of the [Recipes](../recipes.md): droplets that move with and through the gas. The keys are on the [key page](../cases/flow-keys.md) (`spray.*`).

| # | Recipe |
|---|---|
| 242 | [Spray droplets into a flow](#spray-duct) |
| 265 | [Measure the liquid and the droplet sizes at datums](#spray-datums) |
| 280 | [Let the spray push the gas (two-way coupling)](#spray-two-way) |
| 283 | [Cool a gas with droplets (two-way heat)](#spray-cooling) |
| 285 | [Let the turbulence spread the droplets](#spray-dispersion) |
| 292 | [Solve a spray in a steady flow (steady spray mode)](#spray-steady) |
| 306 | [Let the droplets evaporate](#spray-evaporation) |
| 318 | [Break a fuel jet up into drops (KH-ACT, KH-RT)](#spray-breakup) |
| 314 | [Evaporate fuel droplets at high ambient pressure](#spray-a-surface) |
| 321 | [Evaporate a fuel: take a liquid of the library](#spray-liquid-library) |
| 325 | [Start a spray from measured droplets](#spray-brno) |
| 320 | [Check the thermodynamics of a high-pressure spray before the CFD](#spray-a-mixing-line) |
| 319 | [Report the liquid length and the vapour penetration (ECN)](#spray-ecn) |
| 324 | [Give the droplet film real-fluid properties at high pressure](#spray-a-film) |
| 326 | [Let distorted drops drag and evaporate as spheroids](#spray-a-distortion) |
| 323 | [Model a multi-hole injector on a sector (ECN Spray G)](#spray-g-sector) |
| 333 | [Let the drops damp the gas turbulence (two-way turbulence)](#spray-turbulence) |
| 335 | [Let the drops collide](#spray-collisions) |
| 341 | [Describe a multi-hole injector by its holes](#spray-device) |
| 343 | [Show your own injector STL](#spray-device-surface) |
| 338 | [Turn the grid against a jet](#spray-g-grid) |
| 344 | [Pick the time step of a spray that atomizes at the nozzle](#spray-child-size) |
| 363 | [Let colliding drops bounce (Post outcomes)](#spray-collision-post) |
| 369 | [Let the drops meet in collision cells of their own](#spray-collision-cells) |

---

<a id="spray-duct"></a>
## 242. Spray droplets into a flow

**Add a `[spray]` table and an injector file.** Sprays run in transient flow cases (`[time]`). The case names a file with the injectors and the liquid, so the spray can be changed without touching the flow (D246):

```toml
[spray]
include = "spray_duct_injectors.toml"
wall_mode = "stick"          # or "rebound", "escape", "mirror"
write_every = 50             # steps between parcel files
```

```toml
# spray_duct_injectors.toml
[liquid]
name = "water"
density = 998.0

[[injector]]
name = "nozzle"
position = [0.1, 0.395, 0.05]   # inside the mesh (m)
direction = [0.3, -1.0, 0.0]
speed = 10.0                    # m/s
half_angle = 10.0               # degrees; 0 = a line
mass_flow = 0.01                # kg/s
diameter = 150e-6               # m, one size in this version
parcels_per_second = 2000.0     # how finely the flow is sampled
```

Gravity comes from `[gravity]`. `[spray.boundary]` sets the behaviour per boundary patch (default: the `wall_mode` for walls, `escape` for outlets and inlets).

**Run and look.**

```bash
foehn flow examples/cases/flow/spray_duct.toml --out output/spray_duct
python3 -m tools.plot.spray_fields OUT_DIR output/spray_duct   # the picture of gallery 242
```

Open `output/spray_duct/spray/parcels_*.vtkhdf` in ParaView and colour the points by `speed` or `diameter`. `spray.csv` holds the mass book per step; the run stops with SPR005 if it does not close. Check the file only with `foehn flow FILE --check`: SPR001 to SPR003 name a bad injector, a position outside the mesh or an unsupported setting (spray in a case with rotating zones or periodic boundaries).

**Choose the parcels.** More parcels per second give a smoother cloud and cost time in proportion; the droplet physics does not change. A restart continues the parcels bit for bit. See the [gallery entry](../gallery/sprays.md#spray-duct) for what this version leaves out.

---

<a id="spray-datums"></a>
## 265. Measure the liquid and the droplet sizes at datums

**Give the droplets a size distribution.** In the injector file, `size` replaces `diameter` (D242). The distribution is one of the mass over the diameter, cut to `min` and `max`:

```toml
[[injector]]
name = "nozzle"
# position, direction, speed, half_angle, mass_flow, parcels_per_second as in 242
size = { distribution = "rosin-rammler", mean = 120e-6, spread = 2.5, min = 10e-6, max = 400e-6 }
# or: size = { distribution = "log-normal", median = 80e-6, sigma = 1.8, min = 10e-6, max = 300e-6 }
```

`mean` is the Rosin-Rammler size X (63.2 % of the mass in smaller droplets) and `spread` its exponent n; `median` is the mass median diameter and `sigma` the geometric standard deviation of the log-normal form.

**Add datums.** A `[[datum]]` is a named place in the case that the reports use for the gas and the spray (D262). Surfaces count what passes through them, positive along their direction; volumes count what is inside:

```toml
[[datum]]
name = "x06"                 # column prefix in datums.csv
type = "plane"
point = [0.6, 0.2, 0.05]
normal = [1.0, 0.0, 0.0]     # the positive direction
# radius = 0.1               # a disc, or up/width/height for a rectangle
# snap = "faces"             # move it to the nearest mesh faces: exact gas flow

[[datum]]
name = "x16"
type = "surface"
file = "spray_datums_section.stl"   # ASCII or binary STL; normals by the corners

[[datum]]
name = "outlet"
type = "boundary"
boundary = "xmax"

[[datum]]
name = "nozzle_box"
type = "box"                 # also "cylinder" (start, end, radius) and "zone"
lower = [0.0, 0.2, 0.0]
upper = [0.3, 0.4, 0.1]
histogram = { bins = 40, min = 1e-6, max = 1e-3, scale = "log" }
```

`type = "interface"` with `zones = ["a", "b"]` is the face set between two cell zones. A cut plane is exact on the parcels; its gas flow uses the cell values and gradients, and `cut_error` says how far it is from the same plane moved to the faces. Choose `snap = "faces"` when the gas flow must close with the boundaries.

**Run and look.**

```bash
foehn flow examples/cases/flow/spray_datums.toml --out output/spray_datums
python3 -m tools.plot.spray_datums OUT_DIR output/spray_datums   # the picture of gallery 265
```

`datums.csv` gets a row every `[output] monitor_every` steps: per datum the gas columns and the cumulative liquid each way (`liquid_kg_pos`, `liquid_kg_neg`), the net rate of the interval, the crossings, D10, D32 and Dv50. `spray_hist_<datum>.csv` has the droplets and the mass per size bin of each interval. `result.json` has the totals of the run per datum (with D20, D30, D43, Dv10, Dv90, the span, the velocity statistics, a Rosin-Rammler fit, the mean rate and its relative standard error), and the spray in the flow; `spray.csv` the penetration of each injector (`[spray] penetration_fraction = 0.95`, `penetration = "axial"` or `"radial"`). `series.json` lists the field and parcel files for a viewer.

**Ask after the run, without a rerun.** The same numbers are tools (D267), so a script or an agent can ask them:

```bash
echo '{"case": "examples/cases/flow/spray_datums.toml", "restart": "output/spray_datums/restart.h5"}' \
  | foehn call post.datums --input -          # gas through the datums (also for steady runs)
echo '{"case": "examples/cases/flow/spray_datums.toml", "file": "output/spray_datums/spray/parcels_000300.vtkhdf", "datum": "nozzle_box"}' \
  | foehn call post.droplets --input -      # diameters and histogram of the parcels in the box
echo '{"case": "examples/cases/flow/spray_datums.toml", "file": "output/spray_datums/spray/parcels_000300.vtkhdf"}' \
  | foehn call post.penetration --input -     # penetration per injector
```

In Python:

```python
import foehncraft as fc
r = fc.Result.load("output/spray_datums")
r.datums("outlet")            # columns of datums.csv for one datum
r.spray()                     # spray.csv
r.spray_histogram("floor")    # spray_hist_floor.csv
r.series()                    # series.json
p = r.parcels(r.series()[-1]["files"]["parcels"])   # points and fields (needs h5py)
```

A datum with an unknown type, a missing key or a boundary the mesh does not have stops the check with DAT001, one that meets no part of the mesh with DAT002. Datums cannot be used with a sliding mesh in this version. See the [gallery entry](../gallery/sprays.md#spray-datums) for what is checked and what is not there yet.

<a id="spray-two-way"></a>
## 280. Let the spray push the gas (two-way coupling)

**Switch the coupling on.** By default the gas does not feel the droplets (`mode = "one-way"`). With two-way coupling the drag on every parcel is a momentum source of the gas in the cells it passes (D247 to D252):

```toml
[spray.coupling]
mode = "two-way"
passes = "converged"         # or a fixed number of passes, e.g. 1
tolerance = 1e-3             # coupling lag at which the passes stop
max_passes = 3               # at most this many passes in one time step
max_reynolds_change = 0.2    # a sub-step is halved while the particle Reynolds number changes more
```

Two-way coupling needs a transient case with PISO and no solid regions in this version (SPR003); with `[energy]` the heat is coupled too ([recipe 283](#spray-cooling)). Each pass tracks all parcels against one frozen gas state; the gas takes the summed sources implicitly (value plus slope, like a porous sink), so even a dense spray at a time step of many droplet response times stays bounded. With `"converged"` the passes run inside the outer iterations of the time step until the coupling lag, the change of the gas velocity in the cells with parcels against the largest slip or gas change, is below `tolerance`. A step that stops at `max_passes` above it warns once with SPR006 and is counted in `spray.csv`. `passes = 1` is the linearised implicit coupling: cheaper, still conservative and bounded.

**Run and look.**

```bash
foehn flow examples/cases/flow/spray_two_way.toml --out output/spray_two_way
python3 -m tools.plot.spray_two_way OUT_DIR output/spray_two_way   # the picture of gallery 280
```

`spray.csv` gets the columns `passes`, `coupling_lag`, `unconverged`, `gas_force_x/y/z` (N, what the parcels hand to the gas), `gas_power_w` and the health numbers of D254: `coupling_max` (the largest β Δt / (ρ V), how strongly the spray holds the gas of one cell), `coupling_share`, `liquid_fraction_max` and `parcels_per_cell_max`. The ledger in `result.json` books the exchange as a pair (`spray drag` on the gas momentum, `from the gas` on the parcel momentum, pairs 2801 to 2803) that cancels to the bit, and the parcel momentum account closes with the weight, the injection, the walls and what left the flow. See the [gallery entry](../gallery/sprays.md#spray-two-way) for what is checked and what is not there yet.

<a id="spray-cooling"></a>
## 283. Cool a gas with droplets (two-way heat)

**Give the liquid and the injector a temperature.** With `[energy]` in the case, a two-way spray exchanges sensible heat with the gas (SP2b): every parcel carries a temperature that relaxes exactly towards the gas of its cell, with the Nusselt number of Ranz and Marshall, and the heat is a source of the gas energy in value-plus-slope form. In the injector file:

```toml
[liquid]
name = "water"
density = 998.0
specific_heat = 4181.0       # J/(kg K), default 4181 (water)
temperature = 293.15         # K, default 293.15; the temperature of the injected liquid

[[injector]]
name = "nozzle"
# ...
temperature = 293.15         # optional: this injector's liquid temperature, else [liquid] temperature
```

The fluid needs `specific_heat` and `conductivity` (MAT007 otherwise), and the case `[energy]` with an inlet temperature. Then switch the coupling on as in [recipe 280](#spray-two-way):

```toml
[spray.coupling]
mode = "two-way"
tolerance = 1e-3
max_passes = 4               # heat converges more slowly per pass than momentum
```

The passes stop when both lags are below `tolerance`: the velocity lag of recipe 280 and the temperature lag (the change of the gas temperature of a spray cell against the frozen state of the pass, over the largest temperature difference between gas and droplets). The temperature lag falls by a factor of about 7 per pass; when `spray.csv` shows many `unconverged` steps (SPR006), raise `max_passes` to 4 or shorten the step. A one-way spray in a case with `[energy]` heats or cools its droplets too, without acting on the gas.

**Run and look.**

```bash
foehn flow examples/cases/flow/spray_cooling.toml --out output/spray_cooling
python3 -m tools.plot.spray_cooling OUT_DIR output/spray_cooling   # the picture of gallery 283
```

`spray.csv` gets `t_mean_k` (the mass mean droplet temperature) and, two-way, `gas_heat_w` (W, the heat the parcels hand to the gas; negative when they cool it) and `coupling_heat_max` (the largest G Δt / (ρ c_p V) of a cell). The ledger books the exchange as the pair 2804: `spray heat` on the thermal energy account of the energy equation and `from the gas` on the account `parcel enthalpy`, which closes with the injected liquid and what left the flow. The parcel files have a `temperature` field. See the [gallery entry](../gallery/sprays.md#spray-cooling) for what is checked and what is not there yet (no evaporation before SP5).

<a id="spray-dispersion"></a>
## 285. Let the turbulence spread the droplets

**Switch the dispersion on.** By default a parcel sees the mean gas velocity. In a case with a turbulence model, the eddy-lifetime model (D243; Gosman and Ioannides 1983) adds the eddies:

```toml
[spray]
include = "my_injectors.toml"
dispersion = "eddy-lifetime"   # default "none"; "random-walk" follows behind the same key
```

Each parcel then sees a Gaussian fluctuation of variance 2k/3 per component, kept for the eddy interaction time (the eddy lifetime, or the time a heavy droplet needs to cross the eddy). It works one way and two way, with any turbulence model of the code (ε directly from the k-ε family, β* k ω from SST); a laminar case is refused (SPR001). The random numbers are keyed by the parcel and its eddy count, so the run gives the same bits on any number of threads and the two-way passes see the same eddies. In a 2D case the eddies stay in the x-y plane.

**Run and look.**

```bash
foehn flow examples/cases/flow/spray_dispersion.toml --out output/spray_dispersion
python3 -m tools.plot.spray_dispersion OUT_DIR output/spray_dispersion   # the picture of gallery 285
```

The parcel files carry the positions to measure the spread; with `[[datum]]` planes across the duct (recipe 265) `datums.csv` and the histograms give the distribution of the liquid without post-processing. The dispersion depends on the k and ε the turbulence model gives at the nozzle: check that the inlet turbulence does not decay before the spray reaches it (the gallery case uses a viscosity ratio of 2000 for that). Near walls the isotropic eddies overestimate the wall-normal fluctuation, so wall deposition is not reliable with this model. See the [gallery entry](../gallery/sprays.md#spray-dispersion) for what is checked.

<a id="spray-steady"></a>
## 292. Solve a spray in a steady flow (steady spray mode)

**Leave out `[time]` and add `[spray.steady]`.** A spray in a steady case
runs as *streams* (D260, SP2s): each stream stands for a share of the
injector's mass flow, is followed through the frozen gas field, and hands
the gas its sources every `every` outer iterations. All keys have defaults:

```toml
[spray]
include = "spray_cooling_injectors.toml"   # the same injector file as the transient case
wall_mode = "stick"

[spray.coupling]
mode = "two-way"             # "one-way": the streams are followed once at the end, the gas is not touched

[spray.steady]
streams = 400                # streams per injector, default 200
every = 10                   # outer iterations between stream passes, default 10
relaxation = 0.5             # S <- S + relaxation (S_new - S), default 0.5
tolerance = 1e-3             # source residual for convergence, default 1e-3
max_time = 5.0               # s, a stream still in the flow after this is "unfinished" (SPR007), default 10
seed = 0                     # other seeds draw other streams, default 0
```

`parcels_per_second` of the injector file is not needed in the steady
mode. The case needs SIMPLE or SIMPLEC for two-way coupling (SPR003);
`[spray.steady]` with `[time]`, or an injector with `start` or `end`, is
refused (SPR002).

**Run and look.**

```bash
foehn flow examples/cases/flow/spray_cooling_steady.toml --out output/spray_cooling_steady
foehn flow examples/cases/flow/spray_cooling.toml --out output/spray_cooling       # optional, the transient run to compare
python3 -m tools.plot.spray_cooling_steady OUT_DIR output/spray_cooling_steady output/spray_cooling   # the picture of gallery 292
```

The run is converged only when the gas residuals are below `[solver]
tolerance` and the source residual (how far the relaxed sources the gas
holds are from those of the last pass) is below `[spray.steady]
tolerance`; a run that stops at `max_iterations` above it warns with
SPR008. `spray.csv` has one row per pass with the source residual
(`source_residual`, `residual_momentum`, `residual_heat`), the liquid per
boundary in kg/s, the residence times, the force and heat on the gas and
the health numbers (`coupling_max` = β/ṁ, `coupling_heat_max` =
G/(ṁ c_p), with ṁ the gas mass flow through the cell, and
`liquid_fraction_max`). `spray/streams.vtkhdf` holds the paths of the
final pass; in ParaView colour the points by `temperature` or `diameter`.
The log and `result.json` (`spray_steady`) give the summary.

**When it does not converge.** A source residual that stalls above the
tolerance usually means the passes are too far apart for the gas to follow
or the relaxation is too large: lower `relaxation` to 0.3, or `every` to 5.
Strong coupling (health numbers above 1) needs a smaller relaxation. To see
the statistical error of the streams, run again with another `seed` or more
`streams`. See the [gallery entry](../gallery/sprays.md#spray-steady) for
what is checked and how it compares with the transient run.

---

<a id="spray-evaporation"></a>
## 306. Let the droplets evaporate

**Name the vapour and pick a model.** The liquid file says which species
of `[species]` the liquid evaporates into; `[spray.evaporation]` switches
the evaporation on. The case needs `[energy]`, `[species]` with that
species, and two-way coupling:

```toml
# spray_evaporation_injectors.toml
[liquid]
name = "water"
density = 998.0
specific_heat = 4181.0
temperature = 293.15
vapour = "H2O"            # a species of [species]
# saturation = "water"    # IAPWS-95 vapour pressure, the only law so far
# latent_heat = 2.4417e6  # J/kg at latent_heat_temperature = 298.15 K
```

```toml
[species]
names = ["H2O", "O2", "N2"]
initial = { mass_fractions = { H2O = 0.005, O2 = 0.2318, N2 = 0.7632 } }
properties = "fluid"      # the heat balance then closes to round-off

[spray.evaporation]
model = "ranz-marshall"   # or "abramzon-sirignano" (film thickening by the blowing)
# max_mass_change = 0.1          # a sub-step ends when d² falls by more than this share
# max_temperature_change = 2.0   # ... or the droplet temperature changes by more (K)
# min_diameter = 1e-6            # smaller droplets evaporate in one sub-step at the end of their life

[spray.coupling]
mode = "two-way"
passes = "converged"
max_passes = 8            # evaporation couples more slowly than heat alone
```

The same keys work in the steady spray mode (`[spray.steady]`,
[recipe 292](#spray-steady)).

**Run and look.**

```bash
foehn flow examples/cases/flow/spray_evaporation.toml --out output/spray_evaporation
python3 -m tools.plot.spray_evaporation OUT_DIR output/spray_evaporation              # the picture of gallery 306
foehn flow examples/cases/flow/spray_saturation_steady.toml --out output/spray_saturation_steady
python3 -m tools.plot.spray_saturation_steady OUT_DIR output/spray_saturation_steady  # the picture of gallery 306b
```

`spray.csv` adds `evaporated_kg` and `evaporated_parcels` (one way and two
way) and, with two-way coupling, `gas_vapour_kg_s` (the vapour the gas
receives), `vapour_lag_kg_s`, `coupling_vapour_max` (the vapour slope over
the gas mass of the cell per step) and the limiter counts `limiter_cells`
and `limiter_kg_s`. The ledger has the vapour pair 2805 (`spray vapour` on
the gas mass, `to the gas` on the account `parcel mass`) and `spray
vapour` in the species summary. `Y_H2O` is a field of `fields.vtkhdf`.

**When the log warns.** SPR006: the passes of a step ended at `max_passes`
above the tolerance; raise `max_passes` or lower the time step. SPR009: the
saturation limiter of D253 held vapour back because a cell evaporated more
than its gas can take; that is a sign of too few passes or too large a
step (transient) or of a steady run that is not converged. A liquid whose
`vapour` is not a species of the case is refused. See the
[gallery entry](../gallery/sprays.md#spray-evaporation) for what is
checked and what is not.

---

<a id="spray-breakup"></a>
## 318. Break a fuel jet up into drops (KH-ACT, KH-RT)

**Describe the nozzle, not the drops.** Give the injector its hole instead
of a drop size and speed: the parcels then leave as blobs of the diameter
√C_a d_hole at U = ṁ/(ρ_l C_a A_hole), and a measured rate of injection
replaces the constant mass flow (D306):

```toml
# spray_breakup_injectors.toml
[liquid]
name = "diesel-like fuel"
density = 830.0
specific_heat = 2000.0
temperature = 300.0
surface_tension = 0.027    # N/m, needed by the breakup and the dynamic drag
viscosity = 2.5e-3         # Pa s
# vapour_pressure = 5e4    # Pa: the cavities of kh-act may collapse (else they only burst)

[[injector]]
name = "hole"
position = [0.0005, 1e-5, 1e-5]
direction = [1.0, 0.0, 0.0]
half_angle = 6.0
hole_diameter = 1.5e-4
area_contraction = 0.86    # C_a: the blob is √C_a d_hole, the jet speed ṁ/(ρ C_a A)
flow_share = 0.25          # a quarter model of the jet takes a quarter of the liquid
rate = [[0.0, 0.0], [5e-5, 4.46e-3], [1.45e-3, 4.46e-3], [1.5e-3, 0.0]]   # [t, kg/s]
parcels_per_second = 1e6
hole_length = 6e-4         # kh-act: the nozzle turbulence from the hole (Huh and Gosman)
discharge_coefficient = 0.8
# inlet_loss = 0.5         # K_c of the inlet, 0.5 for a sharp edge
# nozzle_k = 240.0         # or give k0 and epsilon0 at the exit directly,
# nozzle_epsilon = 6e8     # for example from a nozzle flow simulation
```

**Pick the model.** `kh-act` lets the turbulence and the cavities the
nozzle leaves in the liquid strip the blobs as well as the
Kelvin-Helmholtz waves, whichever is faster (D307); the children have
the Kelvin-Helmholtz size B₀Λ in either case (D348); `kh-rt` is the
aerodynamic model alone. Both shed child parcels off the blobs and break the
drops beyond the breakup length by Rayleigh-Taylor waves:

```toml
[spray.breakup]
model = "kh-act"       # or "kh-rt"
drag = "dynamic"       # TAB distortion raises the drag towards a disc's
# b0 = 0.61            # KH child radius B0 Λ
# b1 = 10.0            # KH time constant: the one a validation calibrates
# c_rt = 1.0           # RT drop size C_RT Λ_RT
# c_tau = 1.0          # RT time C_τ / Ω_RT
# breakup_length = 5.0 # C_bl: RT beyond C_bl √(ρl/ρg) d_hole (default b1/2)
# shed_fraction = 0.03 # a blob sheds a child per 3 % of its liquid
# c_tcav = 1.0         # kh-act: the rate constant of turbulence and cavitation
```

**Let the gas spread the jet.** In a still chamber use k-ε: the SST model
does not build up the turbulence of the thin gas jet during the injection,
and the spray then runs far ahead (D308 notes). Two-way coupling and the
eddy-lifetime dispersion carry the drops with the gas.

```toml
[turbulence]
model = "k-epsilon"
wall_treatment = "scalable"

[spray]
include = "spray_breakup_injectors.toml"
dispersion = "eddy-lifetime"
penetration_fraction = 0.97    # the penetration in spray.csv: 97 % of the liquid

[spray.coupling]
mode = "two-way"
passes = "converged"
max_passes = 3
```

**Run and look.**

```bash
foehn flow examples/cases/flow/spray_breakup.toml --out output/spray_breakup
foehn flow examples/cases/flow/spray_breakup.toml --out output/spray_breakup_khrt --set spray.breakup.model=kh-rt
python3 -m tools.plot.spray_breakup OUT_DIR output/spray_breakup output/spray_breakup_khrt   # the picture of gallery 318
```

`spray.csv` adds `child_parcels` (children the blobs shed), `rt_breakups`
and, with `kh-act`, `act_children` (the children the turbulence or the
cavitation stripped); `hole.penetration_m` is the distance that holds
`penetration_fraction` of the liquid. The parcel files have `distortion`
with the dynamic drag. A child's `parcel_id` holds its blob's injector and
injection index. Breakup is refused in the steady spray mode. See the
[gallery entry](../gallery/sprays.md#spray-breakup) for what is checked and
what is not.
<a id="spray-a-surface"></a>
## 314. Evaporate fuel droplets at high ambient pressure

**Above about 1 MPa, take the Peng-Robinson surface.** Raoult's law gives
too little fuel vapour at the droplet surface at diesel pressures (1.5 to
1.9 times less for n-dodecane in nitrogen at 6 MPa). Add to the
evaporation:

```toml
[spray.evaporation]
model = "abramzon-sirignano"
equilibrium = "peng-robinson"   # default "raoult"

[spray.evaporation.peng_robinson]
# fuel = "n-dodecane"           # default: [spray.liquid] name
ambient = { N2 = 1.0 }          # far-field mole fractions of N2, O2, CO2, H2O, AR
# k_ij = { N2 = 0.205 }         # per gas; built-in for n-dodecane + N2, else 0
# table_step = 1.0              # K between the nodes of the surface table
```

For a liquid without built-in data give its critical point, acentric
factor and molar mass (all four, then with the κ(ω) of Peng and Robinson,
which can be 10 to 20 % off on the vapour pressure at low temperature):

```toml
[spray.evaporation.peng_robinson]
critical_temperature = 647.096  # K
critical_pressure = 22.064e6    # Pa
acentric_factor = 0.3443
molar_mass = 0.018015           # kg/mol
```

The surface is tabulated once per run at p₀ (open domains; a closed volume
is refused with evaporation). Droplets hotter than the end of the table
(the critical mixing temperature above the fuel's critical pressure) take
no latent heat and are counted in the tally (`critical_parcels`,
`critical_mass`). See the [gallery entry](../gallery/sprays.md#spray-a-surface)
for what is checked and what is not.

---

<a id="spray-a-mixing-line"></a>
## 320. Check the thermodynamics of a high-pressure spray before the CFD

**Run the mixing line of your fuel and gas.** For n-dodecane into
nitrogen (ECN Spray A) the example writes the line and the liquid lengths:

```bash
cargo run --release -p foehncraft-particles --example spray_a_mixing_line -- output/spray_a_mixing_line
python3 -m tools.plot.spray_a_mixing_line output/spray_a_mixing_line output/spray_a_mixing_line
```

`line.csv` has the mixing temperature, the liquid mass fraction, the
density and the energy-balance error over the mixture fraction for the
frozen, the equilibrium and the ideal line; `liquid_length.csv` the
saturation mixture fraction Z_sat of each model and Siebers' liquid length
from it at the measured conditions.

**Your own streams.** Build the mixture as for the droplet surface
(recipe 314), give the ideal-gas enthalpy of every component (the species
data of the flow, read with `foehncraft_particles::ideal_gas_species`) and
the two streams:

```rust
let n2 = ideal_gas_species(&std::fs::read_to_string("crates/foehncraft-flow/data/species-gri30.toml")?, "N2")?;
let fuel = N_DODECANE.ideal_gas;
let h_ig = |i: usize, t: f64| R * t * if i == 0 { fuel.h_over_rt(t) } else { n2.h_over_rt(t) };
let line = Streams { mix: &mix, ideal: &h_ig, p: 6e6, t_fuel: 363.0, t_gas: 900.0, gas };
let z_sat = line.saturation(Mixing::Equilibrium);
```

Compare Z_sat of `Mixing::Equilibrium` and `Mixing::Ideal`: when they
differ by more than a few per cent, the droplet surface needs
`equilibrium = "peng-robinson"`. A liquid length from
`siebers_liquid_length` is a scaling estimate (C_L = 0.62); quote it next
to the measured value, not instead of it. See the
[gallery entry](../gallery/sprays.md#spray-a-mixing-line) for what is
checked and what is not.
<a id="spray-ecn"></a>
## 319. Report the liquid length and the vapour penetration (ECN)

**Add `[spray.ecn]`.** The Engine Combustion Network compares spray
simulations by two numbers defined in its modelling standard
(ecn.sandia.gov, Modeling Standards): the liquid length, the farthest
axial position where the projected liquid volume PLV = ∫ LVF dy across the
spray falls to a threshold, and the vapour penetration, the farthest axial
position where the fuel mass fraction is 0.1 %. With the table, both are
written for every injector and step:

```toml
[spray.ecn]
# plv = [2e-7, 2e-6]          # m: 0.2e-3 and 2e-3 mm³/mm², one liquid length each
# bin = 2.5e-4                # m: square bins along the axis and across
# line_of_sight = [0, 0, 1]   # default: the coordinate axis most nearly across the injector
# vapour_fraction = 1e-3      # fuel mass fraction of the vapour penetration
```

`spray.csv` gets `<injector>.liquid_length_<k>_m` (k counts the `plv`
thresholds from 1) and, with `[spray.evaporation]`,
`<injector>.vapour_penetration_m`; `result.json` lists the same per
injector (`liquid_length`, `vapour_penetration`). The PLV of a bin is the
liquid volume of the injector's parcels in the bin, summed along the line
of sight, over the bin area; the profile along the axis takes the largest
bin across, and the length is interpolated between bin centres. The sums
run in the order of the parcel identities, so the lengths are the same
bits on any number of threads.

The defaults are those of the standard for a round jet such as Spray A.
For a planar (2D) case the PLV is per metre of depth and the thresholds
are yours to choose. With few parcels a single parcel can exceed a
threshold on its own and the length jumps with the farthest parcels: the
standard assumes enough parcels that a bin holds many (check that one
parcel's volume over `bin`² is well below the smaller threshold). A
background of the fuel vapour in the gas counts for the vapour penetration
(the fraction is absolute). Not with `[spray.steady]`. See the
[gallery entry](../gallery/sprays.md#spray-ecn).
<a id="spray-liquid-library"></a>
## 321. Evaporate a fuel: take a liquid of the library

**Name the liquid instead of typing its properties.** In the injector
file:

```toml
[liquid]
library = "n-heptane"          # "water", "n-heptane", "n-dodecane", "iso-octane"
temperature = 300.0            # K, the injected liquid
vapour = "NC7H16"              # its species in [species]: NC7H16, IC8H18, H2O
vapour = "NC7H16"              # its species in [species]: NC7H16, IC8H18, NC12H26, H2O
# property_temperature = 330.0 # K of the constant density and c_l (default: temperature)
```

The library sets the vapour pressure law, the density, the specific heat
and the latent heat; leave `density`, `specific_heat`, `latent_heat` and
`saturation` out (a value given beside `library` is refused). The droplet
heats up to its wet-bulb temperature, so for long-lived droplets set
`property_temperature` near it (n-heptane in hot nitrogen: 330 to 340 K).

**Put the vapour in the gas.** The vapour must be a transported species of
the case, not the closure species:

```toml
[species]
names = ["NC7H16", "N2"]
initial = { mass_fractions = { NC7H16 = 0.0, N2 = 1.0 } }
```

**The film.** `[spray.evaporation] film = "one-third"` is the default: the
gas properties around each droplet are those of the vapour and carrier
mixture at the 1/3 reference state, from the species data. The carrier is
the start composition without the vapour. `film = "far-field"` gives the
first version's behaviour (the gas properties of the cell), which
overstates the evaporation of fuel droplets by 50 to 100 % against
measurements. Case: `examples/cases/flow/spray_heptane_evaporation.toml`;
[gallery entry](../gallery/sprays.md#spray-liquid-library) for what is
checked and what is not.

<a id="spray-brno"></a>
## 325. Start a spray from measured droplets

**When the droplets of a plane are known**, from a phase Doppler or
shadowgraph measurement or from a finer simulation of the atomizer, start
the parcels from them instead of a cone and a size law (D322). Write a
table of one droplet per row, SI units:

```text
# r (m), d (m), axial and radial velocity (m/s); optional u_tangential, weight
r, d, u_axial, u_radial
2.98e-03, 1.164e-05, 7.25, -6.67
5.20e-03, 1.866e-05, 60.0, 13.4
```

and name it in the injector (in the included file, the path is relative to
that file):

```toml
[[injector]]
name = "plane_z20"
position = [0.0, 0.0, 0.0201]   # on the axis, in the measured plane
direction = [0.0, 0.0, 1.0]     # the axis of the spray
mass_flow = 3.5556e-4           # kg/s through the whole plane
samples = "spray_brno_samples.csv"
sector = 90.0                   # a quarter model between symmetry planes
flow_share = 0.25               # = sector / 360
```

Each parcel takes the diameter and velocity of one row, drawn by mass
(`weight`, or equal shares), at the row's radius on a random azimuth, so
the measured link between size, speed and place is kept. Leave out
`speed`, `half_angle`, `diameter`, `size` and `hole_diameter`. With a
`sector` the azimuth runs from the global axis least aligned with
`direction` (x for a spray along z) towards the next one (y): put the
symmetry planes at x = 0 and y = 0. A plane of measured points is
turned into such a table by drawing droplets by their liquid volume and
the ring area of their point; `examples/studies/brno_spray/pda_to_csv.py`
does it for the Brno data. Case: `examples/cases/flow/spray_brno.toml`;
[gallery entry](../gallery/sprays.md#spray-brno) for what is checked and
what is not.


---

<a id="spray-a-film"></a>
## 324. Give the droplet film real-fluid properties at high pressure

**Nothing to do with the Peng-Robinson surface.** With
`equilibrium = "peng-robinson"` and the 1/3-rule film (both D321 and D313
defaults), the film takes the Peng-Robinson density and c_p, Chung's
dense-fluid viscosity and conductivity, and the Riazi-Whitson diffusivity:

```toml
[spray.evaporation]
model = "abramzon-sirignano"
equilibrium = "peng-robinson"
# film = "one-third"            # default; the real-fluid film needs it

[spray.evaporation.peng_robinson]
ambient = { N2 = 1.0 }          # the film's gases
# film = "real-fluid"           # default; "ideal-gas" for the D321 film
# critical_volume = 7.52e-4     # m³/mol, for a liquid given by its critical point
```

At Spray A conditions (6 MPa, 900 K) the real fluid changes the film by 1
to 7 % and the droplet life by about −3 %. At 15 MPa the change is up to
14 %. Compare both with `film = "ideal-gas"` when a result depends on it.

For a liquid given by `critical_temperature`, `critical_pressure`,
`acentric_factor` and `molar_mass`, `critical_volume` sets the molar
volume of Chung's method. The default is Z_c R T_c/p_c with Z_c = 0.291 −
0.080 ω (Pitzer et al. 1955): 2 % below the measured value for n-dodecane,
15 % above it for water. The film properties of one droplet
are written by the example:

```bash
cargo run --release -p foehncraft-particles --example spray_a_film -- output/spray_a_film
python3 -m tools.plot.spray_a_film output/spray_a_film output/spray_a_film
```

See the [gallery entry](../gallery/sprays.md#spray-a-film) for what is
checked and what is not.

---

<a id="spray-a-distortion"></a>
## 326. Let distorted drops drag and evaporate as spheroids

**One key with breakup.** The TAB distortion y of each drop (`[spray.breakup]`,
D303) flattens it (y > 0) or stretches it (y < 0). With the Corrected
Distortion model of Nguyen, Dahms, Pickett and Tagliante (2022) the drop
drags and evaporates as that spheroid, with the internal circulation of a
viscous drop:

```toml
[spray.liquid]
library = "n-dodecane"            # the viscosity μ_l(T) of the drop

[spray.breakup]
model = "kh-act"
drag = "corrected-distortion"     # or "dynamic" (Liu et al. 1993), "sphere"

[spray.evaporation]
model = "abramzon-sirignano"      # or "ranz-marshall"; both take the factor
```

- **Drag:** the viscous sphere of Feng and Michaelides (2001a) for
  λ = μ_l/μ_g, times the shape ratio of Richter and Nikrityuk (2012) in
  l/d = 1 − y/2. A fully flattened drop drags 1.0 to 1.9 times a sphere
  (Re 20 to 500, λ 2.5 to 40), not the 3.6 times of the dynamic drag.
- **Evaporation and heat:** the Sherwood and Nusselt numbers of the sphere
  laws times the viscous-sphere factor (Feng and Michaelides 2001b), the
  shape ratio and the surface ratio S/(π d²). A flattened drop evaporates
  2 to 5 times faster than a solid sphere; a prolate one about as fast as
  a round drop of the same liquid (the internal circulation alone adds 1
  to 46 %).
- **λ** is the liquid's viscosity at the drop temperature (the liquid
  library, else `[spray.liquid] viscosity`) over the gas (film) viscosity.
  The liquid needs a viscosity and a surface tension, as for breakup.

The factors of one drop over y, Re and λ are written by the example:

```bash
cargo run --release -p foehncraft-particles --example spray_a_distortion -- output/spray_a_distortion
```

The Spray A case `examples/cases/flow/spray_a_distortion.toml` runs the
model with `--set spray.breakup.drag=dynamic` as the comparison. See the
[gallery entry](../gallery/sprays.md#spray-a-distortion) for what is
checked and what is not.
<a id="spray-g-sector"></a>
## 323. Model a multi-hole injector on a sector (ECN Spray G)

**Mesh one plume's share.** The plumes of an injector of n equal holes
repeat every 360°/n about its axis, and each plume is mirror-symmetric
about its own plane. A sector of 180°/n with the plume on one side plane
holds half a plume, and both side planes are mirror planes. `wedge` with
`dimension = 3` revolves a box of hexahedra into such a sector (D309): x is
the injector axis, y the radius, and the z intervals become the angular
cells.

```toml
[mesh.generate]
type = "hex"
dimension = 3
lower = [0.0, 3.0e-4, 0.0]   # a slip core of 0.3 mm on the axis
upper = [0.06, 0.05, 1.0]    # 60 mm along the axis, 50 mm radius
cells = [72, 65, 8]          # 8 angular cells
wedge = 22.5                 # 8 holes: 180°/8
segments.x = [ { to = 0.03, cells = 50 }, { to = 0.06, cells = 22, ratio = 3.0 } ]
segments.y = [
  { to = 0.0012, cells = 3, ratio = 2.0 },   # small cells next to the core: no skewed faces
  { to = 0.03, cells = 48 },
  { to = 0.05, cells = 14, ratio = 3.0 },
]

[boundary.ymin]
type = "symmetry"    # the core
[boundary.zmin]
type = "symmetry"    # the plane of the plume
[boundary.zmax]
type = "symmetry"    # half way to the next plume
```

Keep the thin core. Without it the cells next to the axis are prisms whose
faces have a skewness of 1/(4Δφ), 5 with 8 angular cells, and in our trials
they were the first to go unstable once the plumes drew gas to the axis.
By the symmetry of the plumes the gas on the axis moves only along it, and
the slip core imposes just that.

**Put the hole on the plume plane.** The angle runs from −wedge/2 at `zmin`
to +wedge/2 at `zmax`. Put the hole exit a few µm inside the plane
φ = −wedge/2, point the plume along that plane, and give the sector half of
the hole's liquid:

```toml
[[injector]]
name = "hole1"
position = [1.0e-4, 9.81175461e-04, -1.93128751e-04]   # r = 1 mm on φ = −11.25°, 2 µm inside
direction = [0.83867057, 0.53417395, -0.10625380]      # 33° from the axis, in that plane
half_angle = 12.5
flow_share = 0.5           # half of the plume is in the sector
hole_diameter = 1.65e-4
area_contraction = 0.65
rate = [ [0.0, 0.0], [2.7e-5, 1.95e-3], [6.5e-4, 1.72e-3], [7.87e-4, 0.0] ]   # [t, kg/s] of one hole
```

**Measure as the ECN measures.** `tools/plot/spray_g_metrics.py` turns the
sector's parcels into the whole injector (each parcel mirrored about its
plume plane and turned to every plume). `liquid_penetration(x, d, n,
threshold)` projects the liquid volume along the line of sight
perpendicular to the plane of plumes 1 and 5 and returns the farthest
distance along the injector axis above the threshold (0.2e-3 or
2e-3 mm³/mm²). `vapour_penetration(centres, y)` returns the farthest cell
above 0.1 % vapour.

```sh
foehn flow examples/cases/flow/spray_g_g1.toml --out output/spray_g_g1
foehn flow examples/cases/flow/spray_g_g1.toml --set spray.breakup.model=kh-rt --out output/spray_g_g1_khrt
python3 -m tools.plot.spray_g_g1 OUT_DIR output/spray_g_g1 output/spray_g_g1_khrt   # gallery 323
```

Case: `examples/cases/flow/spray_g_g1.toml` with its injector file;
[gallery entry](../gallery/sprays.md#spray-g) for the comparison with the
measurements.

---

<a id="spray-turbulence"></a>
## 333. Let the drops damp the gas turbulence (two-way turbulence)

**Switch it on.** With two-way coupling and the eddy-lifetime dispersion,
the drops can give back to the turbulence what its eddies spend on them
(D330, the term of KIVA-II):

```toml
[spray]
dispersion = "eddy-lifetime"        # needed: the work is that of these eddies

[spray.coupling]
mode = "two-way"                    # needed
turbulence = "dispersion-work"      # default "none"
```

Every sub-step of a parcel books the work of its eddy on the parcel's
response to it, −m (1 − e^{−Δt/τ}) u′·(u′ − v′), which is never positive.
k loses it per unit mass, ε gets c_s (ε/k) of it (c_s = 1.5) and ω of SST
(c_s − 1)(ω/k), all implicit where they take away, so k stays positive.
The steady spray mode refuses the key for now.

**Read what the drops took.** `spray.csv` gets the column
`gas_turbulence_w` (W in each step, negative) and `result.json` the total
`gas_turbulence_work` (J) of the run. `held_eddies` in the tally of
`result.json` counts eddies a parcel kept to the end of a step because it
had already entered 32 in it (`MAX_EDDIES_PER_STEP`): that happens only
where k is nearly gone, as in a spray that has damped its gas turbulence.

**Expect damping, and a lot of it in a dense spray.** Drops much smaller
than the eddies take energy from them at a rate of about 2φ/τ times k, φ
the mass loading of the liquid and τ the response time of the drops. Near
the nozzle of a fuel spray φ is 1 or more, and k decays there within a
fraction of a millisecond, before the gas jet builds its turbulence up;
a two-equation model cannot grow turbulence back from nearly zero k.
Only damping is modelled: the wakes of
drops larger than the eddies, which can add turbulence, are not. On ECN
Spray G (G1, measured before the collision model) this term made the
spray run further ahead of the measurement, so the case leaves it off; see the
[gallery entry](../gallery/sprays.md#spray-turbulence) for the numbers.

```sh
foehn flow examples/cases/flow/spray_g_g1.toml --out output/spray_g_g1
foehn flow examples/cases/flow/spray_g_g1.toml --set spray.coupling.turbulence=dispersion-work \
    --out output/spray_g_g1_work
python3 -m tools.plot.spray_turbulence_two_way OUT_DIR output/spray_g_g1 output/spray_g_g1_work   # gallery 333
```

---

<a id="spray-collisions"></a>
## 335. Let the drops collide

Dense sprays near the nozzle hold many drops per cell, and their
collisions merge drops and slow the core of the spray. Switch them on in
any transient spray:

```toml
[spray.liquid]
surface_tension = 0.02      # N/m; a library liquid brings its own

[spray.collision]
model = "ntc"               # no-time-counter sampling, O'Rourke's outcomes
# coalescence = false       # every collision a grazing separation (a check)
```

The drops of two parcels in one cell meet at the rate q_i q_j σ v Δt / V
(q the drops of a parcel, σ = π (r_i + r_j)², v their relative speed, V the
cell volume). Below O'Rourke's critical impact parameter they coalesce,
above it they separate with their sizes and pull their velocities
together. Mass, momentum and liquid enthalpy stay exact. The injected blobs
of the breakup model do not collide.

Watch `spray.csv`: `collisions`, `coalescences`, `drop_collisions`,
`merged_parcels`, and `collision_capped_cells`, the cell steps whose
candidate pairs reached the cap of 20000 per parcel. A few capped cells
near the nozzle are harmless (in G1, raising the cap from 1000 changed the
penetration by less than 1 %); many mean the cells are too large for the
spray.

```sh
foehn flow examples/cases/flow/spray_g_g1.toml --out output/spray_g_g1
foehn flow examples/cases/flow/spray_g_g1.toml --set spray.collision.model=none --out output/spray_g_g1_nocol
python3 -m tools.plot.spray_g_g1 OUT_DIR output/spray_g_g1 output/spray_g_g1_nocol \
    --label "without collisions" --name spray_g_g1_collisions.png   # gallery 335
```

Case: `examples/cases/flow/spray_g_g1.toml`;
[gallery entry](../gallery/sprays.md#spray-collisions).

---

<a id="spray-device"></a>
## 341. Describe a multi-hole injector by its holes

Give the injector as it is mounted and drilled, and let the solver work out
where each plume starts and where it points. In the injector file (or in
`[spray]` of the case):

```toml
[[device]]
name = "gdi"
mount = [0.0, 0.0, 0.0]        # the tip (m)
axis = [1.0, 0.0, 0.0]         # into the chamber
reference = [0.0, 1.0, 0.0]    # zero of the azimuth (default: the global
                               # axis least aligned with `axis`)

[[device.pattern]]             # a ring of equal holes
count = 8
polar = 33.0                   # degrees from the axis
radius = 1.0e-3                # exit circle (m)
setback = 1.0e-4               # exits this far into the chamber (m)
first_azimuth = -11.25

[[device.hole]]                # or holes one by one
name = "centre"
polar = 0.0

[[injector]]
name = "g"
device = "gdi"
hole = "all"                   # one injector per hole: g/1 ... g/8, g/centre
# speed, rate, hole_diameter, ... as for any injector (per hole)
```

- `hole = "<name>"` takes one hole (a sector model); `offset` (m) moves
  its start, for example a few µm into the domain off a symmetry plane.
- Keep `position` and `direction` beside `device` and `hole` to check a
  hand-written injector against the drawing: SPR010 warns when they
  disagree.
- `clocking` turns the whole ring about the axis, so one file serves
  several mountings.

Check the result before meshing:

```sh
foehn call spray.devices --input - <<< '{"file": "examples/cases/flow/spray_multihole_device_injectors.toml"}'
foehn flow examples/cases/flow/spray_multihole_device.toml --out output/spray_multihole_device
python3 -m tools.plot.spray_multihole_device OUT_DIR output/spray_multihole_device   # gallery 341
```

Case: `examples/cases/flow/spray_multihole_device.toml`;
[gallery entry](../gallery/sprays.md#spray-device).

---

<a id="spray-device-surface"></a>
## 343. Show your own injector STL

Give the device a surface, and every run draws it with the spray
(`spray/devices.vtu`, listed in `series.json` under `static`):

```toml
[[device]]
name = "gdi"
mount = [0.0, 0.0, 0.0]
axis = [1.0, 0.0, 0.0]
geometry = "standard-gasoline"   # or "standard-diesel"; no holes given:
                                 # the tip brings its own (8 at 37°, 7 at 77.5°)
```

Your own tip:

```toml
geometry = "stl"
file = "my_tip.stl"        # relative to this file; ASCII or binary
scale = 1e-3               # the file is in millimetres
# the frame of the file, when it is not tip at the origin, axis +z into
# the chamber, azimuth zero along +x (in the file's own coordinates):
stl_origin = [0.0, 0.0, 12.5]
stl_axis = [0.0, 0.0, -1.0]
stl_reference = [1.0, 0.0, 0.0]
glyph_length = 2e-3        # length of the hole lines (m)
```

- The drawn bores are those of the file; the holes the spray uses are the
  device's `[[device.hole]]` or `[[device.pattern]]` (or the default ones
  of a standard tip; a file of your own needs its holes given). Draw both
  from the same numbers.
- The surface is checked (turned-over, degenerate and duplicate triangles:
  SPR010, drawn as it is); an open surface is fine, it is a picture and
  the flow does not see it (D344).

Look at it before a run:

```sh
foehn call spray.place --input - <<< '{"file": "examples/cases/flow/spray_g_device_injectors.toml", "out": "output/spray_g_device"}'
foehn flow examples/cases/flow/spray_g_g1.toml --set time.end_time=6e-4 --out output/spray_g_g1
python3 -m tools.plot.spray_g_device OUT_DIR output/spray_g_g1   # gallery 343
```

Open `output/spray_g_device/spray/devices.vtu` in ParaView and colour it
by `part` (0 surface, 1 hole axis, 2 cone). Files: `spray_g_device_injectors.toml`,
`assets/injectors/` (README with frame, units, sources);
[gallery entry](../gallery/sprays.md#spray-device-surface).
<a id="spray-g-grid"></a>
## 338. Turn the grid against a jet

A jet that runs along the cell lines of the grid is carried farther than
one that crosses them. To measure it for your case, run the case twice on
a box mesh whose lines you turn with `grid_angle`:

```toml
[mesh.generate]
type = "hex"
dimension = 3
lower = [0.0, 3.0e-4, 0.0]
upper = [0.06, 0.05, 1.0]
cells = [100, 83, 8]
wedge = 22.5          # a sector about the x axis, as for Spray G
grid_angle = 33.0     # the cell lines at 33° in the (x, r) plane; -12 turns them 45° away
```

The cells are uniform with the spacings `(upper − lower)/cells`, cut by
the sides of the box; a cut piece of less than 0.3 of a cell joins its
neighbour and becomes a polyhedron. It needs hexahedra, no grading,
segments, `refine`, `perturb`, `split`, `remove` or `zones`. The injector
cannot be turned inside a sector (its axis is the sector's axis and its
side planes are the mirror planes of the plumes), so turn the grid, not
the injector. Choose the angle from the plume direction: the plume of
Spray G leaves at 33° from the axis, so 33° puts one family of lines along
it and −12° puts both 45° across it.

```sh
foehn flow examples/cases/flow/spray_g_g1_grid_aligned.toml --out output/spray_g_grid_aligned
foehn flow examples/cases/flow/spray_g_g1_grid_turned.toml --out output/spray_g_grid_turned
foehn flow examples/cases/flow/spray_g_g1.toml --out output/spray_g_g1        # the third curve
python3 -m tools.plot.spray_g_grid OUT_DIR output/spray_g_grid_aligned output/spray_g_grid_turned output/spray_g_g1
```

Cases: `examples/cases/flow/spray_g_g1_grid_aligned.toml` and
`spray_g_g1_grid_turned.toml`; [gallery entry](../gallery/sprays.md#spray-g-grid).

<a id="spray-child-size"></a>
## 344. Pick the time step of a spray that atomizes at the nozzle

With `model = "kh-act"` the nozzle turbulence can strip a blob within a few
microseconds into children of the Kelvin-Helmholtz size B₀Λ (D348), drops
of a few µm right at the nozzle. Such drops follow the gas within a
microsecond, and the passes of the two-way coupling then need a small time
step: G1 of Spray G diverges near 0.1 ms at 2 µs and runs at 1 µs; Spray A
(150 MPa, 22.8 kg/m³) needs 0.5 µs, the step of published set-ups.

```toml
[time]
scheme = "bdf2"
time_step = 1e-6     # Spray G G1; 5e-7 for Spray A

[spray]
write_every = 50     # keep the parcel output at 50 µs when the step halves
```

Once the injector stops, no new blobs atomize at the nozzle and the step
can grow again. G1 steps up after the end of injection (0.787 ms) and
saves a restart there, so the late spray can be rerun on its own:

```toml
[run]
restart_times = [9e-4]

[[schedule]]
at_time = 9e-4
set = { time.time_step = 2e-6 }
```

```sh
foehn flow examples/cases/flow/spray_g_g1.toml --out output/g1_late \
  --restart output/spray_g_g1/restart_000900.h5
```

On G1 the 2 µs leg from 0.9 to 1.2 ms matches the run at 1 µs throughout
(same penetrations, liquid mass within 0.2 %) and takes about 2 minutes.

Signs that the step is too large, before the run diverges: in the log the
residuals of the passes alternate between two values instead of falling,
and the Courant number climbs step by step (from 0.3 to 1.7 in 20 steps on
G1 at 2 µs). Check the first 0.1 ms with `--set time.end_time=1e-4` before
a long run.

```sh
foehn flow examples/cases/flow/spray_g_g1.toml --out output/spray_g_g1
python3 -m tools.plot.spray_g_g1 OUT_DIR output/spray_g_g1
```

Case: `examples/cases/flow/spray_g_g1.toml`;
[gallery entry](../gallery/sprays.md#spray-child-size).

---

<a id="spray-collision-post"></a>
## 363. Let colliding drops bounce (Post outcomes)

In a dense spray O'Rourke's outcomes (D331) merge most slow collisions.
The outcomes of Post and Abraham (D378) let drops bounce at small Weber
numbers and separate reflexively in fast head-on collisions:

```toml
[spray.collision]
model = "ntc"
outcomes = "post"    # default "orourke"
```

Use it for sprays into a gas at elevated pressure with many fine drops,
as in Spray G (6 bar): it brought the drops 15 mm from the tip from 14.5
to 13.0 µm and the liquid lead at 0.8 ms from 9 to 4 %. The bounces and
reflexive separations are counted in the step tally; `spray.csv` keeps
its columns. Compare with the run without collisions
(`--set spray.collision.model=none`) to see how much the collisions do.

```sh
foehn flow examples/cases/flow/spray_g_g1.toml --out output/spray_g_g1 --end-time 9e-4
python3 -m tools.plot.spray_g_g1 OUT_DIR output/spray_g_g1 output/spray_g_g1_orourke \
    --label "O'Rourke outcomes" --name spray_collision_post.png
```

Case: `examples/cases/flow/spray_g_g1.toml`;
[gallery entry](../gallery/sprays.md#spray-collision-post).

---

<a id="spray-collision-cells"></a>
## 369. Let the drops meet in collision cells of their own

By default drops collide when they share a gas cell, so the collisions
depend on the mesh. Give the collisions a lattice of their own, turned
and shifted at random every step (D379):

```toml
[spray.collision]
model = "ntc"
outcomes = "post"
cell = 6e-4    # cube edge in m; leave out to use the gas cells
```

Start with a cube about the size of the gas cells in the dense part of
the spray. Smaller cubes merge more drops; compare two sizes before you
trust a drop size (in Spray G 0.6 mm gave 10.8 µm 15 mm from the tip,
0.2 mm gave 14.0 µm). The run stays bitwise reproducible on any thread
count.

```sh
foehn flow examples/cases/flow/spray_g_g1.toml --out output/spray_g_g1 --end-time 9e-4
# the comparison: a copy of the case without the cell line
foehn flow spray_g_g1_gas.toml --out output/spray_g_g1_gas --end-time 9e-4
python3 -m tools.plot.spray_g_g1 OUT_DIR output/spray_g_g1 output/spray_g_g1_gas \
    --label "flow cells" --name spray_collision_cells.png
```

Case: `examples/cases/flow/spray_g_g1.toml`;
[gallery entry](../gallery/sprays.md#spray-collision-cells).
