# Porous regions

Part of the [Recipes](../recipes.md): flow through porous matrices. The keys are on the [key page](../cases/flow-keys.md) (`regions.<name>.porous.*`).

| # | Recipe |
|---|---|
| 273 | [Put a porous block or a honeycomb into a flow](#porous-flow) |
| 275 | [Give the matrix a temperature of its own: a packed bed or a heat store](#porous-two-temperatures) |

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
porosity = 0.4           # optional here (needed with the two temperatures, 275)
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
support (scalars, turbulence, rotating zones, the GPU, axes turned against
the mesh).

---

<a id="porous-two-temperatures"></a>
## 275. Give the matrix a temperature of its own: a packed bed or a heat store

**Name the material of the matrix.** With `[energy]`, a porous region needs
`solid`, a material with density, specific heat and conductivity; the
region then carries a gas and a matrix temperature in every cell (D271).

```toml
[materials.ceramic]
kind = "solid"
density = 2300.0
specific_heat = 900.0
conductivity = 2.0

[regions.bed]
kind = "porous"
zones = ["bed"]
[regions.bed.porous]
viscous = 1.0e4
porosity = 0.4                      # needed: ε of the gas, 1 − ε of the matrix
solid = "ceramic"
heat_transfer_coefficient = 10.0    # h (W/(m² K)) between gas and matrix
specific_area = 60.3                # a (1/m), wetted surface per volume
solid_conductivity = 0.0            # optional: effective k of the matrix, one value or three along the axes
solid_initial = 300.0               # optional: start of the matrix (number or expression); default the gas start
source = "chosen values"
```

A `honeycomb` gives all of it from the geometry: h = Nu k_gas / d_h
(`nusselt`, default 2.98, laminar square ducts), a = 4 ε / d_h, and the
conductivity of the brick walls, (1 − ε) k_s along the channels and
(1 − √ε) k_s across (D272):

```toml
[regions.stack.porous]
honeycomb = { channel = 6.9e-3, pitch = 13.6e-3, length = 0.15 }
solid = "ceramic"
```

The gas in the region takes ε ρ c and ε k of `[fluid]`; the velocity
stays the superficial one. The matrix is adiabatic at the faces of the
region (contact with a housing is PM3).

**Run and look.**

```bash
foehn flow examples/cases/flow/porous_schumann.toml --out output/schumann
foehn flow examples/cases/flow/porous_storage_charge.toml --out output/storage
python3 -m tools.plot.porous_energy_fields OUT_DIR output/schumann output/storage
```

The field `matrix_temperature` is in `fields.vtkhdf` (the gas temperature
outside porous regions), and the restart file carries it. The heat balance
in the log and in `heat.csv` books the matrix as an account of its own
(`<region> matrix` in the storage). POR001 names a missing `solid`,
`porosity`, `heat_transfer_coefficient` or `specific_area`; POR002 refuses
`[energy.cycle]`, `[energy.acceleration]` and a density that follows the
temperature in porous cases.
