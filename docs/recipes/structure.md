# Structures

Part of the [Recipes](../recipes.md): the linear structural solver `foehncraft-structure` (L53, stage MS5). A structural case is its own TOML file ([keys](../cases/structure.md#structural-case-file-keys)), run through the tool registry.

| # | Recipe |
|---|---|
| 295 | [Run a structural case of beams and clamp springs, and fit the clamp to a measured frequency](#structure-case) |
| 297 | [Model a flat reed as a shell clamped along an edge](#structure-plate) |
| 299 | [Run a structure in time, with damping, load histories and temperatures](#structure-transient) |

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
## 299. Run a structure in time, with damping, load histories and temperatures

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


---

<a id="structure-solid"></a>
## 300. Model a thick part as a solid, from Gmsh or as a block

1. Mesh it in Gmsh with second-order elements: ten-node tetrahedra, or
   twenty-node hexahedra with `Mesh.SecondOrderIncomplete = 1` (the
   27-node ones are refused). Give the volume a physical group, and every
   surface that is held or loaded a physical group of its own. Save as
   `.msh` 4.1, ASCII or binary:

   ```sh
   gmsh stop.geo -3 -order 2 -format msh41 -o reed_stop.msh
   ```

   (`python3 -m tools.mesh.reed_stop_msh` writes the stop of the gallery
   without Gmsh.)

2. Add it as a part. `mesh` is relative to the case file; `volumes`
   picks physical volumes (default: all):

   ```toml
   [[solids]]
   name = "stop"
   mesh = "reed_stop.msh"
   volumes = ["stop"]
   material = "steel"
   ```

   Or generate a block, with the sides `xmin`, `xmax`, `ymin`, `ymax`,
   `zmin`, `zmax` as its surfaces:

   ```toml
   [[solids]]
   name = "seat"
   block = { origin = [0.0, 0.0, 0.0], size = [0.03, 0.01, 0.004], divisions = [12, 4, 2], element = "hex20" }
   material = "steel"
   ```

3. Hold and load it on its named surfaces. A pressure on a solid pushes
   into it (a fluid's pressure); `force_per_area` is a traction in global
   axes:

   ```toml
   [[supports]]
   part = "stop"
   surface = "mount"
   fix = ["all"]

   [[pressures]]
   part = "stop"
   surface = "underside"
   value = 3.0e5
   ```

   `[[temperatures]]` takes a uniform `change` for a solid. Points,
   modes and transients work as for beams and plates.

4. Read the results: `result.json` lists the solid parts with their
   surfaces and gives `static.max_von_mises` (the stress at the nodes,
   averaged over the elements around each); `structure.vtk` is an
   unstructured grid of quadratic cells with `displacement`, `von_mises`
   and the mode shapes.


---

<a id="structure-contact"></a>
## 302. Let a reed touch its seat and its stop

1. Give every seat and stop a `[[contacts]]` plane: a `point`, the
   `normal` pointing to the side the reed stays on, the `part` whose
   nodes may touch it, and `offset` = half the thickness of a plate or a
   beam, so that its face touches, not its mid-surface:

   ```toml
   [[contacts]]
   name = "seat"
   part = "reed"
   point = [0.0, 0.0, -0.1e-3]
   normal = [0.0, 0.0, 1.0]
   offset = 0.1e-3

   [[contacts]]
   name = "stop"
   part = "reed"
   point = [0.0, 0.0, 0.3e-3]
   normal = [0.03331483, 0.0, -0.99944491]   # rises 1 mm over 30 mm
   offset = 0.1e-3
   ```

   For a solid, `surface` limits the plane to the nodes of a named
   surface. A node that starts inside a plane is an error.

2. A static case with contacts is solved exactly (to `[contact]
   tolerance`, 10⁻⁹ of the model size). A transient needs `method =
   "direct"`; take a step that resolves the motion between impacts (here
   5 µs, 1/1300 of the first period) and `output_every = 1` to see every
   impact, and leave `rho_infinity` below 1 so that the high frequencies
   an impact excites die out.

3. Run it and read the results: `result.json` gives per plane the force
   along its normal (static: also the nodes in contact), and the largest
   penetration; `transient.csv` has a column `contact_<name>` per plane:

   ```sh
   target/release/foehn call structure.run --input examples/cases/tools/structure_reed_valve_stop.json
   python3 -m tools.plot.structure_contact output/reed_valve_stop reed_valve.png
   ```

   `[contact] penalty` (default 1000) scales the penalty stiffness of a
   transient; a larger one gives less penetration and more active-set
   changes.
