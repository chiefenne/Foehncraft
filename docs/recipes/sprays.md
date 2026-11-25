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
Kelvin-Helmholtz waves, whichever is faster (D307); `kh-rt` is the
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
