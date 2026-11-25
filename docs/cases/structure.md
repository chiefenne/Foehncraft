# Structural cases (`examples/cases/structure/`)

Part of the [Case catalogue](../cases.md): structural cases of the linear finite element solver `foehncraft-structure` (L53, stage MS5). A structural case is a TOML file of its own, run by the tool `structure.run` (`foehn call structure.run --input <file>.json`, inputs in [`examples/cases/tools/`](tools.md)); it writes `result.json` and `structure.vtk` into its output folder.

| Case | What it is | Mesh | Discussed in |
|---|---|---|---|
| `reed_beam.toml` | a reed valve of spring steel, 30 × 6 × 0.2 mm, as 40 Timoshenko beams on a rigid clamp: deflection under 2 kPa, eight modes | built from `[[beams]]` | [gallery 295](../gallery/structure.md#reed-beam), [recipe 295](../recipes/structure.md#structure-case) |
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
| `[materials.<name>]` | `young` (Pa), `poisson`, `density` (kg/m³): isotropic linear elastic |
| `[sections.<name>] shape` | `rectangle` (`width` along local y, `thickness` along local z), `circle` (`diameter`) or `general` (`area`, `iy`, `iz`, `torsion`) |
| `[sections.<name>] kappa_y`, `kappa_z` | shear coefficients (default Cowper 1966: 10(1+ν)/(12+11ν) for rectangles, 6(1+ν)/(7+6ν) for circles, 5/6 for `general`) |
| `[[beams]]` | a straight line `from` → `to` cut into `elements` equal beams with `section`, `material`, `theory` (`timoshenko`, default, or `euler-bernoulli`) and `y_axis` (direction of the section's local y; default: local z as close to global z as possible); lines that share an end point share its node; `name` for line loads |
| `[[supports]]` | `at` a node, `fix` a list of `ux`, `uy`, `uz`, `rx`, `ry`, `rz` or `all` |
| `[[clamps]]` | the mounting of a reed (D290) at a node: `normal` (axis of the lift) and `line` (axis of the clamp line); every other unknown of the node is held; `k_translation` (N/m per m) along the normal and `k_rotation` (N·m/rad per m) about the line, times `width` (default: the rectangle's width); a spring left out is rigid |
| `[[springs]]` | `at` a node, `translation` and `rotation` stiffnesses to ground along the global axes (N/m, N·m/rad) |
| `[[masses]]` | `at` a node, `mass` (kg): a translational point mass |
| `[[loads]]` | `at` a node, `force` (N) and `moment` (N·m), global axes |
| `[[line_loads]]` | `beam` (a line's name), `force_per_length` (N/m, global axes), uniform |
| `[points]` | `name = [x, y, z]`: nodes whose displacements (static) and mode shapes are reported |
| `[output] vtk` | write `structure.vtk` (default true): legacy VTK lines with `displacement`, `rotation` and `mode_<n>`, `mode_<n>_rotation` scaled to a largest translation of 1 |

`result.json` (`foehncraft.structure-result`, version 1) holds the sizes,
the mass, every clamp with its stiffness or `"rigid": true`, the static
result (largest displacement, strain energy, reaction force, the named
points) and the modal result (frequency, ω, effective mass fractions in
x, y, z and the residual of each mode; the Sturm check; the mass-normalised
shapes at the named points), and the warnings (`STR003`).
