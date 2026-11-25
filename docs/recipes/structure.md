# Structures

Part of the [Recipes](../recipes.md): the linear structural solver `foehncraft-structure` (L53, stage MS5). A structural case is its own TOML file ([keys](../cases/structure.md#structural-case-file-keys)), run through the tool registry.

| # | Recipe |
|---|---|
| 295 | [Run a structural case of beams and clamp springs, and fit the clamp to a measured frequency](#structure-case) |
| 297 | [Model a flat reed as a shell clamped along an edge](#structure-plate) |
| MMM | [Run a structure in time, with damping, load histories and temperatures](#structure-transient) |

---

<a id="structure-case"></a>
## 295. Run a structural case of beams and clamp springs, and fit the clamp to a measured frequency

1. Write the case: materials, a section, the beam lines, the clamp and
   what to compute.

   ```toml
   [analysis]
   static = true
   modes = 8

   [materials.spring_steel]
   young = 206e9
   poisson = 0.3
   density = 7850.0

   [sections.reed]
   shape = "rectangle"
   width = 6e-3
   thickness = 0.2e-3

   [[beams]]
   name = "reed"
   from = [0.0, 0.0, 0.0]
   to = [0.03, 0.0, 0.0]
   elements = 40
   section = "reed"
   material = "spring_steel"

   [[clamps]]
   at = [0.0, 0.0, 0.0]
   normal = "z"          # the reed lifts along z
   line = "y"            # the clamp line runs across it
   k_rotation = 40.0     # N·m/rad per metre of clamp line; leave out for rigid

   [[line_loads]]
   beam = "reed"
   force_per_length = [0.0, 0.0, 12.0]

   [points]
   tip = [0.03, 0.0, 0.0]
   ```

2. Run it with a JSON input that names the case and the output folder:

   ```sh
   echo '{"case": "examples/cases/structure/reed_beam_clamp_springs.toml", "out": "output/reed"}' > reed.json
   target/release/foehn call structure.run --input reed.json
   ```

   The frequencies go to standard error, the result envelope with
   `result.json` to standard output. `output/reed/structure.vtk` opens in
   ParaView: colour by `displacement` or warp by `mode_1`.

3. Check what the run took as given: `result.json` lists every clamp with
   its total stiffness or `"rigid": true`, and `modal.sturm.ok` is `true`
   when the factorisation's Sturm count agrees with the modes found (a
   `STR003` warning otherwise: raise `modes`).

4. Fit the clamp to a measured first natural frequency (from a tap test or
   the supplier's data sheet):

   ```sh
   echo '{"case": "examples/cases/structure/reed_beam_clamp_springs.toml", "spring": "rotation", "mode": 1, "frequency": 150.0}' > fit.json
   target/release/foehn call structure.fitclamp --input fit.json
   ```

   The answer gives `k_per_length` to put into `k_rotation`, and
   `frequency_rigid`, the limit of a rigid clamp. `STR004` means the
   measurement lies above that limit: check the free length, the
   thickness and the modulus before blaming the clamp.

5. For a free part (no supports), give `[analysis] shift` a negative
   value of about −0.01 (2π f₁)² so that K − σM can be factorised; the six
   rigid-body modes come out at 0 Hz first.

---

<a id="structure-plate"></a>
## 297. Model a flat reed as a shell clamped along an edge

1. Give the plate as a parallelogram: a corner, two edges and the
   divisions along them. The normal is `edge_a × edge_b`; a positive
   pressure pushes along it.

   ```toml
   [[plates]]
   name = "reed"
   origin = [0.0, -6e-3, 0.0]
   edge_a = [0.03, 0.0, 0.0]     # along the reed
   edge_b = [0.0, 12e-3, 0.0]    # across it: normal +z
   divisions = [40, 16]
   thickness = 0.2e-3
   material = "spring_steel"
   element = "mitc4"             # or "mitc3+"
   ```

2. Clamp it along the edge: `along` takes the two ends of a segment and
   holds every node on it. The springs are per metre of clamp line and
   are shared by the length each node stands for; the clamp line's length
   is the segment's, so leave out `width`.

   ```toml
   [[clamps]]
   along = [[0.0, -6e-3, 0.0], [0.0, 6e-3, 0.0]]
   normal = "z"
   line = "y"
   k_rotation = 40.0
   k_translation = 2.0e7
   ```

   `[[supports]]` take `along` in the same way, for example a simply
   supported edge with `fix = ["uz"]`.

3. Load it with a pressure (Pa, along the normal) or a force per area
   (N/m², global axes):

   ```toml
   [[pressures]]
   plate = "reed"
   value = 2000.0
   ```

4. Run it with `structure.run` as in recipe 295. `structure.vtk` holds the
   shell as polygons: in ParaView, *Warp By Vector* with `mode_2` shows
   the twist. Use 2 to 4 elements across a reed of a few millimetres for
   the bending modes and 8 to 16 for the twist and the modes across the
   width.

---

<a id="structure-transient"></a>
## MMM. Run a structure in time, with damping, load histories and temperatures

1. Add `[transient]` to a structural case. The run starts from rest at
   t = 0; `time_step` sets the accuracy (20 to 50 steps per period of the
   highest mode that matters):

   ```toml
   [transient]
   end_time = 0.02
   time_step = 2.5e-5
   method = "direct"          # generalised-α; or "modal"
   rho_infinity = 0.9         # 1: no numerical damping; lower damps the highest modes
   snapshots = [0.001, 0.002] # fields into structure.vtk
   ```

   `method = "modal"` with `modes = n` superposes the lowest n modes,
   each integrated exactly: fast, and right when the load does not excite
   the modes left out.

2. Damp it. `ratios` gives a damping ratio at two frequencies and the
   run finds C = αM + βK; `rayleigh = [α, β]` gives the factors; the
   modal method also takes `modal_ratio` for every mode:

   ```toml
   [damping]
   ratios = [[154.0, 0.02], [1002.0, 0.02]]
   ```

3. Give loads a history. A table is linear between its points and holds
   its last value; a sine takes `amplitude`, `frequency` (Hz), `phase`
   (degrees) and `offset`. Loads without `history` are constant.

   ```toml
   [histories.pulse]
   table = [[0.0, 0.0], [0.0005, 1.0], [0.001, 0.0]]

   [[pressures]]
   plate = "reed"
   value = 3000.0
   history = "pulse"
   ```

4. Load it with temperatures. The material needs `expansion` (1/K);
   `change` is the rise from the stress-free state, `difference` the top
   minus bottom across the thickness (plates: top along the normal; beams:
   +z minus −z, and `difference_y` across the width):

   ```toml
   [materials.spring_steel]
   young = 206e9
   poisson = 0.3
   density = 7850.0
   expansion = 1.2e-5

   [[temperatures]]
   part = "reed"
   change = 80.0
   difference = 6.0
   ```

   Temperatures work in static and transient runs, and take a `history`.

5. Read the results: `result.json` has `transient` (times, the named
   points' displacements, kinetic and strain energy, the largest
   displacement and when), `transient.csv` the same as columns.

