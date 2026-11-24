# Porous regions

Part of the [Recipes](../recipes.md): flow through porous matrices. The keys are on the [key page](../cases/flow-keys.md) (`regions.<name>.porous.*`).

| # | Recipe |
|---|---|
| 273 | [Put a porous block or a honeycomb into a flow](#porous-flow) |

---

<a id="porous-flow"></a>
## 273. Put a porous block or a honeycomb into a flow

**Mark the zone as porous and give its matrix.** A porous region is a fluid
region with a resistance; the velocity is the superficial one.

```toml
[regions.block]
kind = "porous"
zones = ["block"]
[regions.block.porous]
viscous = 1.0e6          # K⁻¹ (1/m²): one value, or three along the axes
inertial = 10.0          # F (1/m); the sink is (μ K⁻¹ + ½ ρ |u| F) u
porosity = 0.4           # optional here (used by the two temperatures, PM2)
source = "Ergun fit of the bed"
```

Direction-dependent values need the axes of the region, the first two
principal directions in the mesh frame (for now along the mesh axes):

```toml
[regions.brick]
kind = "porous"
zones = ["brick"]
axes = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]   # channels along x
[regions.brick.porous]
honeycomb = { channel = 6.9e-3, pitch = 13.6e-3, length = 0.15 }
```

`honeycomb` replaces `viscous`, `inertial` and `porosity`: ε = (channel /
pitch)², K⁻¹ along = (f_D Re / 2) / (ε d_h²) with f_D Re = 56.9, the
entrance loss F = K(∞) / (ε² L) with K(∞) = 1.43 when `length` is given,
and `blocking` (default 10⁴) times K⁻¹ across.

**Run and look.**

```bash
foehn flow examples/cases/flow/porous_darcy_block.toml --out output/porous_block
foehn flow examples/cases/flow/porous_honeycomb_channel.toml --out output/along
foehn flow examples/cases/flow/porous_honeycomb_channel.toml --set regions.brick.porous.honeycomb.blocking=1.0 --out output/iso
foehn flow examples/cases/flow/porous_darcy_block.toml --set 'mesh.generate.cells=[50, 4, 1]' --out output/block50
foehn flow examples/cases/flow/porous_darcy_block.toml --set 'mesh.generate.cells=[200, 4, 1]' --out output/block200
python3 -m tools.plot.porous_fields OUT_DIR output/along output/iso output/block50 output/porous_block output/block200
```

The inlet pressure is in `result.json` (`solves[0].boundaries`, `mean_pressure`
of `xmin`). Check the file only with `foehn flow FILE --check`: POR001 names
a missing or contradicting matrix, POR002 a setting this stage does not
support (energy, scalars, turbulence, rotating zones, the GPU, axes turned
against the mesh).
