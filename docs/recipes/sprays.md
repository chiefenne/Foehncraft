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
| 319 | [Report the liquid length and the vapour penetration (ECN)](#spray-ecn) |

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
