# Lagrangian sprays

Part of the [Recipes](../recipes.md): droplets that move with and through the gas. The keys are on the [key page](../cases/flow-keys.md) (`spray.*`).

| # | Recipe |
|---|---|
| 242 | [Spray droplets into a flow](#spray-duct) |
| 265 | [Measure the liquid and the droplet sizes at datums](#spray-datums) |
| 280 | [Let the spray push the gas (two-way coupling)](#spray-two-way) |

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

Two-way coupling needs a transient case with PISO and no energy equation in this version (SPR003). Each pass tracks all parcels against one frozen gas state; the gas takes the summed sources implicitly (value plus slope, like a porous sink), so even a dense spray at a time step of many droplet response times stays bounded. With `"converged"` the passes run inside the outer iterations of the time step until the coupling lag, the change of the gas velocity in the cells with parcels against the largest slip or gas change, is below `tolerance`. A step that stops at `max_passes` above it warns once with SPR006 and is counted in `spray.csv`. `passes = 1` is the linearised implicit coupling: cheaper, still conservative and bounded.

**Run and look.**

```bash
foehn flow examples/cases/flow/spray_two_way.toml --out output/spray_two_way
python3 -m tools.plot.spray_two_way OUT_DIR output/spray_two_way   # the picture of gallery 280
```

`spray.csv` gets the columns `passes`, `coupling_lag`, `unconverged`, `gas_force_x/y/z` (N, what the parcels hand to the gas), `gas_power_w` and the health numbers of D254: `coupling_max` (the largest β Δt / (ρ V), how strongly the spray holds the gas of one cell), `coupling_share`, `liquid_fraction_max` and `parcels_per_cell_max`. The ledger in `result.json` books the exchange as a pair (`spray drag` on the gas momentum, `from the gas` on the parcel momentum, pairs 2801 to 2803) that cancels to the bit, and the parcel momentum account closes with the weight, the injection, the walls and what left the flow. See the [gallery entry](../gallery/sprays.md#spray-two-way) for what is checked and what is not there yet.
