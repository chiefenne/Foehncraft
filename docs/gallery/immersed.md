# Immersed bodies

Part of the [Example gallery](../gallery.md): bodies in the flow without a body-fitted mesh (L13; [staged plan](../requirements/l13-immersed-boundary-plan.md), D225 to D234).

| # | Example | Stage | Feature |
|---|---|---|---|
| 270 | [A cylinder put into a box mesh by Brinkman penalisation](#penalised-cylinder) | L13 stage 1 | `[immersed]`, `[[immersed.body]]`, solid fraction, body force, `solid-fraction` field |

---

<a id="penalised-cylinder"></a>
## 270. A cylinder put into a box mesh by Brinkman penalisation

<img src="../figures/gallery/immersed_cylinder_re20.png" alt="Four panels of the cylinder in a channel at Re 20. Top left: the speed on a plain box mesh with the cylinder as a penalised body, the exact circle dashed. Top right: the speed of the same case on the body-fitted mesh; the two fields look alike. Bottom left: the vorticity of the penalised run, red below and blue above the cylinder. Bottom right: a close-up of the solid fraction of each cell around the cylinder, 1 inside, 0 outside and partial values in a ring of cells along the dashed circle." width="800" height="346">

`immersed_cylinder_re20.toml` is the steady benchmark 2D-1 of Schäfer &
Turek (1996), the case of [gallery 15](flow.md#cylinder), on a plain box
mesh of 440 × 82 hexahedra (cell width h = D/20) that knows nothing of the
cylinder. The cylinder is one `[[immersed.body]]`: every cell gets its solid
fraction χ (the part of its volume inside the body, from 64 sample points per
tetrahedron of the cell), and the momentum rows get the Brinkman drag
−μ χ u / K with the permeability K = 10⁻⁸ m² (Angot, Bruneau & Fabrie 1999).
The drag is implicit and sits inside the momentum row before the
under-relaxation, so the Rhie–Chow flux sees it and the pressure does not
decouple in the solid. The force on the body is the momentum the sink takes
out of the fluid, Σ μ χ V u / K.

| Mesh | Cells | C_D | C_L | Δp (Pa) |
|---|---|---|---|---|
| box, h = D/10 (220 × 41) | 9 020 | 5.421 | 0.0163 | 0.1203 |
| box, h = D/20 (the case file) | 36 080 | 5.581 | 0.0147 | 0.1187 |
| box, h = D/40 (880 × 164) | 144 320 | 5.571 | 0.0140 | 0.1179 |
| body-fitted, [gallery 15](flow.md#cylinder) | 24 576 | 5.581 | 0.0104 | 0.1175 |
| reference range | | 5.57 – 5.59 | 0.0104 – 0.0110 | 0.1172 – 0.1176 |

Drag and pressure difference reach the reference range at h = D/20 to D/40;
the lift stays 27 % above it at D/40 and falls slowly (about first order).
SIMPLEC and the coupled solver give the same values to all printed digits
(D/10: 542 and 75 iterations; D/40: SIMPLEC does not converge in 3 000
iterations, the coupled solver needs 288).

**What this is for, and its flaws (D226).** Penalisation is the cheap first
look and the topology tool of L13, not the accurate immersed wall: the wall
is smeared over about one cell, so wall shear stress and lift are poor, and
the result depends on K and on the mesh width. Plane Poiseuille flow between
two penalised walls that do not lie on cell faces converges at about first
order (mass flow error 20, 8.5 and 2.3 % on 16, 32 and 64 cells across the
channel; with K = 10⁻⁴ m², √K about a third of a cell, the walls leak and
the flow is 3 % too large). Turbulence models, the energy equation and the
other settings that need a sharp wall are refused (IMM102). Sharp walls on
cut cells are stage 2.

Reproduce (about a minute each):

```sh
cargo build --release -p foehncraft-cli
target/release/foehn flow examples/cases/flow/immersed_cylinder_re20.toml --out output/imm
target/release/foehn flow examples/cases/flow/cylinder_re20.toml --out output/fitted
pip install matplotlib h5py numpy
python3 -m tools.plot.immersed_fields docs/figures/gallery output/imm output/fitted
```

Recipe: [270](../recipes/immersed.md#penalised-body). Case file: `examples/cases/flow/immersed_cylinder_re20.toml`.
