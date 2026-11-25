# Structures

Part of the [Example gallery](../gallery.md): the linear structural solver of the crate `foehncraft-structure` (L53, stage MS5 of the [moving solids plan](../assessments/moving-solids-reed-valves-plan.md)), made for reed valves and other parts that the flow moves. Part 1 has Euler–Bernoulli and Timoshenko beams, clamp springs per length of the clamp line (D290) with the fit to a measured frequency, static analysis and modes by shift-invert Lanczos with a Sturm check. Shells, solid elements, dynamics and contact follow in the next parts of MS5; the coupling to the flow is MS6.

| # | Example | Stage | Feature |
|---|---|---|---|
| 295 | [A reed valve as a beam: modes, and the clamp spring against the frequency equation](#reed-beam) | MS5 | `structure.run`, `structure.fitclamp`, `[[beams]]`, `[[clamps]]` |

---

<a id="reed-beam"></a>
## 295. A reed valve as a beam: modes, and the clamp spring against the frequency equation

<img src="../figures/gallery/structure_reed_modes.png" alt="The first four modes of a reed valve 30 mm long, 6 mm wide and 0.2 mm thick, clamped at one end, drawn as the deformed strip and coloured by the displacement: mode 1 at 184 Hz bends the reed once, mode 2 at 1152 Hz has one node along the length, mode 3 at 1745 Hz twists the strip about its axis, mode 4 at 3225 Hz bends with two nodes" width="990" height="726">

**The reed.** `examples/cases/structure/reed_beam.toml` is a reed of
spring steel (E = 206 GPa, ν = 0.3, ρ = 7850 kg/m³), 30 mm free length,
6 mm wide and 0.2 mm thick, as 40 Timoshenko beams along x, lifting along
z. The clamp at x = 0 holds every unknown of the root; with no springs
given it is rigid, and `result.json` says so (`"rigid": true` for the
translation and the rotation). The case asks for the static deflection
under 2 kPa over the width (12 N/m) and the first eight modes:

```sh
cargo build --release -p foehncraft-cli
target/release/foehn call structure.run --input examples/cases/tools/structure_reed_beam.json
```

```text
reed_beam: 8 modes, Sturm check passed
  mode   1      183.887139 Hz
  mode   2     1152.160832 Hz
  mode   3     1745.491401 Hz
  mode   4     3225.002430 Hz
  ...
reed_beam: static, largest displacement 1.474581e-3 m
```

The tip deflects 1.4746 mm, the deflection of the cantilever under a
uniform load qL⁴/(8EI) + qL²/(2κGA) to the printed digits (shear adds
0.005 %). The flapping modes 1, 2 and 4 are those of the clamped–free beam
(β_n L = 1.875, 4.694, 7.855, Blevins 1979). Mode 3 is the first torsion
mode, which the picture shows because every node of the beam carries its
6 mm section, moved by the axis and turned by the rotation. The Sturm
check counts the eigenvalues of K − σM below and above the modes found
from the signs of the pivots of its factorisation: no mode was skipped.

<img src="../figures/gallery/structure_reed_clamp.png" alt="The first two flapping frequencies of the reed against the stiffness of a rotational clamp spring k L over EI from 0.1 to 1000 on log axes: the solver's marks lie on the lines of the frequency equation, mode 1 rises from 28 Hz to the 184 Hz of the rigid clamp, mode 2 from 800 Hz to 1152 Hz" width="814" height="506">

**The clamp spring (D290).** A real reed is held by a clamp plate that
gives a little: `[[clamps]]` takes a translational stiffness along the
lift (`k_translation`, N/m per metre of clamp line) and a rotational one
about the clamp line (`k_rotation`, N·m/rad per metre), times the clamp
width (the 6 mm of the section). The picture sweeps the rotational spring
from soft (kL/EI = 0.1) to stiff (1000): `structure.run` with 40
Euler–Bernoulli beams lies on the frequency equation of the beam with a
rotational root spring (w(0) = 0, EI w''(0) = k w'(0), free tip) within
1.3·10⁻⁷ for both flapping modes, and both approach the rigid clamp from
below. `reed_beam_clamp_springs.toml` mounts the reed on
k_rotation = 40 N·m/rad per m and k_translation = 2·10⁷ N/m per m: the
first mode falls from 184 Hz to 152 Hz.

**The fit.** `structure.fitclamp` finds the spring that gives a measured
frequency. A stiffer spring raises every natural frequency, so the root is
unique; the fit brackets it by factors of ten and refines it by the
Illinois method on ln k:

```sh
target/release/foehn call structure.fitclamp --input examples/cases/tools/structure_fit_clamp.json
```

```text
clamp 0 rotation: k = 3.693688e1 per metre of clamp line (2.216213e-1 in total), mode 1 at 150.000000 Hz (rigid clamp 183.709967 Hz), 11 modal solves
```

A measured frequency at or above the rigidly clamped model cannot be met
by any spring; the fit stops with `STR004` and says the model is too soft
(free length, thickness or modulus).

**Verification** (`crates/foehncraft-structure/tests/beams.rs`):

| Check | Exact solution | Result |
|---|---|---|
| Tip load on a cantilever, 1, 7 and 40 elements | PL³/(3EI) (+ PL/(κGA)), tip slope PL²/(2EI) | within 10⁻⁸ (rounding of a reed whose axial stiffness is 10⁸ times its bending one) |
| Uniform load, every node | qx²(6L² − 4Lx + x²)/(24EI) (+ shear) | within 10⁻¹⁰ of the tip value |
| Cantilever modes 1 to 5, Euler–Bernoulli | β_n L of Blevins (1979) | within 2·10⁻⁵ on 40 elements, fourth order (3.8 to 4.2); first mode 61.3 % effective mass |
| Pinned round bar, L/d = 10, Timoshenko | the frequency equation of the pinned Timoshenko beam (Han, Benaroya & Wei 1999) | within 10⁻⁵ on 120 elements, second order; below Euler–Bernoulli |
| Clamp springs, static | PL³/(3EI) + PL²/k_r + P/k_t | within 10⁻¹⁰ |
| Rotational clamp spring, kL/EI = 0.3 to 300 | the frequency equation solved by bisection | within 10⁻⁸ |
| Clamp fit | the spring of the frequency equation | recovered within 10⁻⁵ in fewer than 40 solves |
| Free–free beam, negative shift | six rigid modes at 0, β₁L = 4.730 | within 10⁻⁶; the rigid modes carry the whole mass |
| Unsupported static case, a support off the mesh | – | `STR002` naming the free unknown, `STR001` |

The Timoshenko element (interdependent interpolation, Friedman & Kosmatka
1993) has the exact static stiffness, but its consistent mass converges
at second order in the frequency of a thick beam; on a thin reed it is the
Euler–Bernoulli element to rounding. The pictures are made with
`python3 -m tools.plot.structure_reed` (the commands are in its header).
