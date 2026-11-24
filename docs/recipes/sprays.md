# Lagrangian sprays

Part of the [Recipes](../recipes.md): droplets that move with and through the gas. The keys are on the [key page](../cases/flow-keys.md) (`spray.*`).

| # | Recipe |
|---|---|
| 242 | [Spray droplets into a flow](#spray-duct) |

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
