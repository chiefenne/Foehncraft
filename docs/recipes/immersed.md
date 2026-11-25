# Immersed bodies

Part of the [Recipes](../recipes.md): put bodies into the flow without meshing them (L13). The keys are on the [key page](../keys.md) (`immersed.*`).

| # | Recipe |
|---|---|
| 270 | [Put a body in the flow without meshing it](#penalised-body) |

---

<a id="penalised-body"></a>
## 270. Put a body in the flow without meshing it

**Add an `[immersed]` table with one `[[immersed.body]]` per body** to a laminar, isothermal flow case on any mesh. The mesh needs no boundary for the body; the body may cross cells and faces anywhere:

```toml
[immersed]
permeability = 1e-8        # K (m²): sink μ χ V / K; no default
# samples = 4              # sample points per direction and tetrahedron (1 to 16)

[[immersed.body]]
name = "cylinder"
shape = "cylinder"         # "cylinder", "sphere" or "box"
centre = [0.2, 0.2, 0.0]   # a point on the axis
radius = 0.05
# axis = [0.0, 0.0, 1.0]   # default z
# length = 0.1             # default endless

[[immersed.body]]
name = "block"
shape = "box"
lower = [0.6, 0.1, -1.0]
upper = [0.7, 0.2, 1.0]

[output]
derived = ["solid-fraction"]   # χ in fields.vtkhdf
```

**Choose K.** Smaller K is closer to a no-slip wall and stiffer. Keep √K well below a cell width: about 10⁻⁴ h² or less (the gallery case uses 4·10⁻⁴ h²). SIMPLEC takes any K; the coupled solver converges from rest with a sink up to about 10⁵ times the viscous diagonal (μ/K below about 10⁵ μ/h²), above that its linear solver can break down. Check the effect of K by halving it once.

**Run and read the forces.**

```sh
target/release/foehn flow examples/cases/flow/immersed_cylinder_re20.toml --out output/imm
```

The log lists every body with its solid volume on the mesh against the exact volume (IMM104; IMM103 warns when the mesh does not resolve the body). `result.json` has, per solve, `immersed`: name, solid volume, exact volume, cells and the force of the fluid on the body (N). `[verification] benchmark = "schafer-turek-2d1"` takes the body named `cylinder` (or `benchmark_boundary`) when the mesh has no boundary of that name.

**Limits.** Laminar and low-Re only, and isothermal: a turbulence model, the energy equation, passive scalars, gravity, regions, rotating zones and sprays are refused (IMM102). The wall is smeared over about a cell: drag and pressure come out well, wall shear stress and lift do not ([gallery 270](../gallery/immersed.md#penalised-cylinder)). Use it for a first look and topology; a sharp wall comes with the cut cells of L13 stage 2.
