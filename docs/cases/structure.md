# Structural cases (`examples/cases/structure/`)

Part of the [Case catalogue](../cases.md): structural cases of the linear finite element solver `foehncraft-structure` (L53, stage MS5). A structural case is a TOML file of its own, run by the tool `structure.run` (`foehn call structure.run --input <file>.json`, inputs in [`examples/cases/tools/`](tools.md)); it writes `result.json` and `structure.vtk` into its output folder.

| Case | What it is | Mesh | Discussed in |
|---|---|---|---|
| `reed_beam.toml` | a reed valve of spring steel, 30 × 6 × 0.2 mm, as 40 Timoshenko beams on a rigid clamp: deflection under 2 kPa, eight modes | built from `[[beams]]` | [gallery 295](../gallery/structure.md#reed-beam), [recipe 295](../recipes/structure.md#structure-case) |
| `reed_plate.toml` | a reed petal of spring steel, 30 × 12 × 0.2 mm, as 40 × 16 MITC4 shells clamped along an edge on clamp springs: lift under 2 kPa, six modes | built from `[[plates]]` | [gallery 297](../gallery/structure.md#reed-plate), [recipe 297](../recipes/structure.md#structure-plate) |
| `reed_plate_pulse.toml` | the reed petal struck by a 1 ms pressure pulse, 20 ms in time with Rayleigh damping, four snapshots | built from `[[plates]]` | [gallery 299](../gallery/structure.md#reed-dynamics), [recipe 299](../recipes/structure.md#structure-transient) |
| `reed_plate_hot.toml` | the reed petal 80 K warmer and 6 K hotter on top: thermal bending towards the seat | built from `[[plates]]` | [gallery 299](../gallery/structure.md#reed-dynamics) |
| `reed_stop.toml` | a curved reed stop, 25 mm on a 40 mm radius, 8 × 1.5 mm, as 160 twenty-node hexahedra read from Gmsh: bolted at one end, 0.3 MPa on its underside, von Mises stress, four modes | `reed_stop.msh` (Gmsh 4.1, second order; `python3 -m tools.mesh.reed_stop_msh`) | [gallery 300](../gallery/structure.md#reed-stop), [recipe 300](../recipes/structure.md#structure-solid) |
| `reed_beam_clamp_springs.toml` | the same reed on rotational and translational clamp springs per length of the clamp line (D290); the input of the clamp fit | built from `[[beams]]` | [gallery 295](../gallery/structure.md#reed-beam), [recipe 295](../recipes/structure.md#structure-case) |

<a id="structural-case-file-keys"></a>
## Structural case-file keys

All values in SI units. Points (`at`, `from`, `to`) are coordinates; a
point that names a node must lie on one (within 10⁻⁶ of the model size),
otherwise `STR001` names the nearest node.

| Key | Meaning |
|---|---|
| `[structure] name` | name of the case (default: the file name) |
| `[analysis] static` | solve K u = f (default false) |
| `[analysis] modes` | number of modes (0: no modal analysis) |
| `[analysis] shift` | shift σ (rad²/s²): the modes closest to it; 0 for the lowest modes of a supported structure, a negative value (about −0.01 (2π f₁)²) for a free one |
| `[materials.<name>]` | `young` (Pa), `poisson`, `density` (kg/m³): isotropic linear elastic; `expansion` (1/K) for `[[temperatures]]` |
| `[sections.<name>] shape` | `rectangle` (`width` along local y, `thickness` along local z), `circle` (`diameter`) or `general` (`area`, `iy`, `iz`, `torsion`) |
| `[sections.<name>] kappa_y`, `kappa_z` | shear coefficients (default Cowper 1966: 10(1+ν)/(12+11ν) for rectangles, 6(1+ν)/(7+6ν) for circles, 5/6 for `general`) |
| `[[beams]]` | a straight line `from` → `to` cut into `elements` equal beams with `section`, `material`, `theory` (`timoshenko`, default, or `euler-bernoulli`) and `y_axis` (direction of the section's local y; default: local z as close to global z as possible); lines that share an end point share its node; `name` for line loads |
| `[[plates]]` | a flat parallelogram `origin` + a·`edge_a` + b·`edge_b` (a, b from 0 to 1) cut into `divisions` = [n_a, n_b] shell elements of `thickness` and `material`; `element` = `mitc4` (default, quadrilaterals) or `mitc3+` (triangles, alternating diagonals); normal `edge_a` × `edge_b`; plates and beams that share a point share its node; `name` for pressures |
| `[[solids]]` | a part of quadratic solid elements (three translations per node) with `name` and `material`: `mesh` = a Gmsh `.msh` 4.1 file of second order relative to the case file (ten-node tetrahedra, twenty-node hexahedra; six-node triangles and eight-node quadrilaterals of a physical surface name a surface) with `volumes` = the physical volumes to take (default all); or `block` = { `origin`, `size`, `divisions` = [n_x, n_y, n_z], `element` = `hex20` (default) or `tet10` (six per cell) } with the surfaces `xmin` … `zmax`; a solid shares no nodes with beams and plates |
| `[[pressures]]` | `plate` (a plate's name), `value` (Pa, along the plate's normal) and `force_per_area` (N/m², global axes), uniform; or `part` and `surface` (a solid's named surface), `value` pushing into the solid and `force_per_area` |
| `[[supports]]` | `at` a node, `along` = [[x, y, z], [x, y, z]] every node on that segment, or `part` and `surface` every node of a solid's named surface; `fix` a list of `ux`, `uy`, `uz`, `rx`, `ry`, `rz` or `all` |
| `[[clamps]]` | the mounting of a reed (D290) `at` a node (a beam) or `along` a segment (a shell edge, the springs shared by the length each node stands for, the clamp line as long as the segment): `normal` (axis of the lift) and `line` (axis of the clamp line); every other unknown of the node is held; `k_translation` (N/m per m) along the normal and `k_rotation` (N·m/rad per m) about the line, times `width` (default: the rectangle's width); a spring left out is rigid |
| `[[springs]]` | `at` a node, `translation` and `rotation` stiffnesses to ground along the global axes (N/m, N·m/rad) |
| `[[masses]]` | `at` a node, `mass` (kg): a translational point mass |
| `[[loads]]` | `at` a node, `force` (N) and `moment` (N·m), global axes |
| `[[line_loads]]` | `beam` (a line's name), `force_per_length` (N/m, global axes), uniform |
| `[[temperatures]]` | `part` (a beam line, a plate or a solid; a solid takes `change` only), `change` (K) from the stress-free state, `difference` (K) across the thickness (plates: top along the normal minus bottom; beams: +z minus −z), `difference_y` (beams, +y minus −y), `history` |
| `history` (in `[[loads]]`, `[[line_loads]]`, `[[pressures]]`, `[[temperatures]]`) | the name of a `[histories.<name>]` that scales the load in time; none: constant |
| `[histories.<name>]` | `table` = [[t, factor], …] (linear, the last value held) or `amplitude`, `frequency` (Hz), `phase` (degrees), `offset`: a sine |
| `[transient]` | `end_time`, `time_step` (s), `method` (`direct`: generalised-α, default; `modal`), `rho_infinity` (0 to 1, default 0.9), `modes` (modal), `output_every` (default: at most 1000 kept steps), `snapshots` (times of whole fields in `structure.vtk`); starts from rest |
| `[damping]` | `rayleigh` = [α, β] (C = αM + βK) or `ratios` = [[f1, ζ1], [f2, ζ2]]; `modal_ratio` (modal method, every mode) |
| `[points]` | `name = [x, y, z]`: nodes whose displacements (static) and mode shapes are reported |
| `[output] vtk` | write `structure.vtk` (default true): legacy VTK lines (beams) and polygons (shells) with `displacement`, `rotation` and `mode_<n>`, `mode_<n>_rotation` scaled to a largest translation of 1; with solids an unstructured grid of quadratic cells (VTK types 24 and 25) that adds `von_mises` |

`result.json` (`foehncraft.structure-result`, version 1) holds the sizes (`elements`, `beam_elements`, `shell_elements`, `solid_elements`), the solid parts (`solids`: element, mesh source, surfaces),
the mass, every clamp with its node count, its stiffness or `"rigid": true`, the static
result (largest displacement, strain energy, reaction force including the forces of the springs to ground, the named
points, and with solids `max_von_mises`) and the transient result (times, the named points' displacements, kinetic and strain energy, the largest displacement and its time, the snapshot times; `transient.csv` has the same as columns), the modal result (frequency, ω, effective mass fractions in
x, y, z and the residual of each mode; the Sturm check; the mass-normalised
shapes at the named points), and the warnings (`STR003`).
