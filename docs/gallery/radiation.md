# Thermal radiation

Part of the [Example gallery](../gallery.md): thermal radiation by Monte Carlo transport of energy bundles (step 29, L15). R1a covers grey absorbing-emitting gas on a prescribed temperature field, grey diffuse walls, symmetry planes and open boundaries. Bundles follow the cell walk with exact path lengths per cell and carry integer tallies, so a call gives bitwise the same result on any number of threads and ranks and on the GPU (D255 to D259, D269). The specification is the [radiation page](../requirements/step29-radiation.md). The coupling to the energy equation and the `[radiation]` keys of the case file come with R1b.

| # | Example | Stage | Feature |
|---|---|---|---|
| 259 | [A box furnace with cooled tubes, and the cases with exact solutions](#furnace) | R1a | `foehncraft_radiation::Domain::solve`, `Settings::points`, `output::write_vtkhdf` |
| 281 | [Paths of energy bundles in the box furnace, and the net flux into every wall](#furnace-paths) | R1a | `foehncraft_radiation::Domain::trace` |

---

<a id="furnace"></a>
## 259. A box furnace with cooled tubes, and the cases with exact solutions

<img src="../figures/gallery/radiation_furnace_plane.png" alt="Mid-depth plane of a box furnace 4 m long and 2.5 m high with seven tubes under the roof. Top: the prescribed flame temperature, which rises from the burner on the left to 1950 K and leans up towards the outlet. Bottom: the net radiative source of the gas, strongly negative (the gas emits more than it absorbs) in the flame core and close to zero elsewhere" width="800" height="937">

**The furnace (D255).** `examples/meshes/furnace.py` builds the chamber
with Gmsh: 4 m long, 2 m deep and 2.5 m high, about 58 500 tetrahedra.
Seven tubes of 0.16 m diameter cross it under the roof, the burner is a
disc of 0.5 m in the front wall, and the flue gas leaves through an
opening under the roof of the back wall. The example `furnace` prescribes
the flame: the temperature peaks at 1 950 K about 1.2 m from the burner
and falls to the 1 250 K of the flue gas, and the absorption coefficient
(soot) rises from 0.3 1/m to 1.1 1/m in the core. The walls are grey
(ε = 0.75) at 1 150 K, the tubes grey (ε = 0.85) at 620 K, the burner
disc grey (ε = 0.8) at 1 100 K, and the opening is an open boundary.

With 16 million bundles the call takes about 5 minutes on 4 CPU threads.
Every bundle carries a whole number of quanta (the emitted power over
2⁵²), so the energy account closes exactly:

| Energy account | kW |
|---|---|
| emitted by the gas (4κσT⁴V) and the walls (εσT⁴A) | 10 243.2 |
| absorbed by the gas | 4 712.2 |
| absorbed by the walls, the tubes and the burner | 5 450.7 |
| escaped through the opening | 80.2 |
| removed by Russian roulette (net) | 0.001 |

| Boundary | Area (m²) | Incident (kW) | Net into the wall (kW) | Net (kW/m²) |
|---|---|---|---|---|
| tubes | 6.94 | 891.5 | 708.3 ± 0.5 | 102.0 |
| walls | 40.25 | 5 582.5 | 1 192.8 ± 1.1 | 29.6 |
| front wall | 4.80 | 640.9 | 123.3 ± 0.4 | 25.7 |
| burner | 0.20 | 31.7 | 12.4 ± 0.1 | 63.5 |
| opening (open) | 0.48 | 80.2 | – | – |

<img src="../figures/gallery/radiation_furnace_walls.png" alt="Top: the incident radiative flux on the side wall of the furnace, highest (about 190 kW per square metre) on the wall opposite the flame core, lower under the roof next to the cold tubes. Bottom: the net flux into three of the tubes against the angle around the tube, highest on the underside that faces the flame" width="800" height="780">

The tubes take 708 kW, about 102 kW/m² on average. The tubes
furthest downstream take the most (from 95 kW/m² on the second tube to 110 kW/m² on the seventh, as the flame leans up towards the outlet), and every tube takes most
on its underside, which faces the flame. The side wall sees the flame core
with up to 190 kW/m² and less under the tubes, which are cold and so send
little back.

**Verification (section 3 of the specification).** Every case runs in
`crates/foehncraft-radiation/tests/cases.rs` against an exact solution.
The statistical error σ is the batch-mean error of the call, which is
conservative because the scrambled Sobol points stratify all the bundles
of an emitter while the batches are only random subsets of them.

| Case | Exact solution | Result |
|---|---|---|
| 1. Isothermal grey slab between cold black plates, τ = 0.01 to 10 | wall flux σT⁴ (1 − 2 E₃(τ)), the source of every layer from E₂ | wall flux within 0.14 % for every τ; every layer source within 3σ |
| 2. Isothermal enclosure, black and grey walls | no net source, no net wall flux | mean \|source\| 1.6 to 2 % of the emission, mean z² 0.7 |
| 3. Transparent gas between opposed black plates | configuration factor of parallel rectangles (Howell's catalogue) | 0.19983 against 0.19982, and 0.50866 against 0.50899 |
| 4. Energy account, threads, ranks, GPU data path | emitted = absorbed + escaped + roulette, to the quantum | closes exactly; bitwise the same on 1, 2 and 4 threads, on 1 to 4 ranks (also with MPI in CI) and on the device path |
| 5. Cube of isothermal grey gas (κL = 1), cold black walls | flux into each wall face, the gas emission integrated over the five other walls (`exact::cube_wall_flux`, adaptive Gauss–Legendre to 1e-10) | 600 wall faces of 10³ hexahedra: mean z² 0.74, the largest \|z\| 3.2; every wall mean within 0.08 % on hexahedra and within 0.17 % on perturbed polyhedra |
| 6. Sobol against pseudo-random emission | cases 1 and 3 | below; at 512 000 bundles Sobol is 13 times (slab) and 8 times (plates) more accurate |

<img src="../figures/gallery/radiation_sobol.png" alt="Root-mean-square error of the slab wall flux and of the plate exchange against the number of bundles, for scrambled Sobol and for pseudo-random emission points on log axes: the Sobol errors are lower and fall faster than the square-root line that the pseudo-random errors follow" width="770" height="506">

**Sobol against pseudo-random emission (case 6, the benchmark rule of
section 4).** The position and the direction of every bundle come from
one Owen-scrambled Sobol point per emitter (Burley 2020); a wall face
picks its triangle with the first coordinate and uses the rest of that
coordinate inside the triangle, so the position and the direction of a
wall bundle form one stratified point in four dimensions.
`Settings::points = Points::Random` gives the plain Monte Carlo estimate
with keyed pseudo-random numbers. Over 16 independent calls per point,
the error of the slab wall flux falls as N^-0.70 with the Sobol points
and as N^-0.52 with pseudo-random ones, and the error of the plate
exchange as N^-0.63 against N^-0.50. At 512 000 bundles the Sobol points
are 13 times (slab) and 8 times (plates) more accurate, which is the
accuracy of 180 and 60 times as many pseudo-random bundles.

**Speed.** On a cloud machine with 4 cores, one million bundles in the
furnace take 84 s on one thread, 44 s on two and 20 s on four, about
84 µs per bundle on one thread. A bundle there crosses the chamber
several times between the grey walls before Russian roulette ends it. On
MPI ranks with one thread each, the polyhedral case of the example
`mpi_radiation` (2 431 cells, one million bundles) takes 22.1 s on one
rank, 11.8 s on two and 6.2 s on four, with 1.8 million hand-overs
between the ranks, and the tallies are bitwise those of one rank. The
time on a GPU comes with the CPU = GPU check of the Colab notebook
`tools/notebooks/radiation_colab.ipynb`.

```
pip install gmsh matplotlib numpy
python3 examples/meshes/furnace.py --size 0.15 --out output/furnace.msh
cargo run --release -p foehncraft-radiation --example furnace -- \
    output/furnace.msh output/furnace 16000000
cargo run --release -p foehncraft-radiation --example sobol_study > output/sobol.csv
python3 -m tools.plot.radiation_furnace docs/figures/gallery output/furnace output/sobol.csv
cargo test --release -p foehncraft-radiation --test cases -- --nocapture
```

`output/furnace/furnace.vtkhdf` opens in ParaView with the source, its
error, the temperature and the absorption coefficient in the cells, and
the incident and net flux with its error on every boundary.

[↑ back to top](#thermal-radiation)

---

<a id="furnace-paths"></a>
## 281. Paths of energy bundles in the box furnace, and the net flux into every wall

<img src="../figures/gallery/radiation_paths_3d.png" alt="The box furnace in 3D, seen from above the near side wall, with the seven cooled tubes under the roof. The floor, the far side wall and the front wall with the burner are coloured by the net radiative flux into the wall, brightest on the floor and the side wall opposite the flame and on the burner disc, dark on the side wall around the tube ends. Sixty bundle paths cross the chamber as straight lines between wall reflections, cyan for bundles from the gas and green for bundles from the walls, each fading as its energy is absorbed" width="800" height="512">

**Paths from the call itself.** The furnace of [gallery 259](#furnace)
with the same mesh, flame and walls and the same call of 16 million
bundles. `Domain::trace` takes the arguments of `Domain::solve` and the
bundles to follow (`pick`: how many of each emitter, counted from its
first). It emits these bundles as the call does, with the same emitters,
weights and Sobol points, and moves them with the same kernel, whose
tracing hook does nothing in a call and on the GPU. A traced path is
therefore the path that bundle takes in the call, and the energy it
leaves in the cells and on the walls is its share of the call's result:
the test `tests/trace.rs` traces every bundle of a call of 20 000 and
finds the absorbed power of every cell, the net flux of every wall face
and the energy account of the call again, to the quantum.

The example `furnace` with `PATHS = 60` draws the emitters in proportion
to their power (systematic sampling over the cells, then the wall faces),
so the 60 paths are a fair sample of the call: 40 start in the gas and 20
on the walls. Eleven of the 40 gas bundles start where the gas is hotter
than 1 600 K, which is 5 % of the volume but 23 % of the emission of the
gas (4κσT⁴ per unit volume). Every point of a path is where the bundle
leaves a cell (155 cells per path on average), is reflected at a grey
wall, or ends.

**A path fades with the energy it still carries.** In every cell the gas
takes the share 1 − exp(−κℓ) of it; a grey wall absorbs ε of what
arrives (0.75 for the refractory, 0.85 for the tubes) and sends the rest
back diffusely. The median bundle reaches its first wall with about half
of its energy and leaves it with an eighth, crosses the chamber three or
four more times and ends after 6.4 m, when Russian roulette takes it at a
thousandth of its energy. One of the 60 leaves through the opening; in
the whole call 0.8 % of the emitted power does.

<img src="../figures/gallery/radiation_paths_plane.png" alt="Side view of the furnace, x from 0 to 4 m and z from 0 to 2.5 m, with the gas temperature on the mid-depth plane in Rainbow Desaturated: the flame leaves the burner on the left, peaks at 1950 K and leans up towards the outlet on the right. The 60 bundle paths are drawn over it as white lines from their start points, bright while they carry most of their energy and faint after the first wall" width="800" height="533">

<img src="../figures/gallery/radiation_paths_walls.png" alt="The net radiative flux into the six walls of the furnace on one Inferno scale, in two rows: floor, roof and front wall with the burner disc; the side wall y = 0, the side wall y = 2 m and the back wall with the opening. The floor and both side walls are brightest opposite the flame core, the roof is dark behind the tubes, the side walls have dark rings around the tube ends, and the burner disc is bright" width="800" height="387">

**The net flux into every wall** from the same call: the floor takes
36 kW/m² on average, the side walls 34, the back wall 36, the front wall
with the burner 27 and the roof 10. The floor and the side walls see the
flame core and take up to about 65 kW/m² opposite it. The roof lies
behind the cold tubes (620 K), which take 102 kW/m² themselves and send
little back; next to them the side walls emit more than they receive
(down to −20 kW/m² in the dark rings around the tube ends). The two side
walls agree within the noise, as the furnace is symmetric about its
mid-depth plane.

The paths do not depend on the number of bundles of the call: bundle k
of an emitter takes the same Sobol point, so a call of 200 000 bundles
traces the same 60 paths, in the same points, in a few seconds. The
call of 16 million bundles takes about 6 minutes on 4 threads, tracing
the 60 bundles 0.01 s.

```
pip install gmsh matplotlib numpy
python3 examples/meshes/furnace.py --size 0.15 --out output/furnace.msh
cargo run --release -p foehncraft-radiation --example furnace -- \
    output/furnace.msh output/furnace 16000000 60
python3 -m tools.plot.radiation_paths docs/figures/gallery output/furnace
cargo test --release -p foehncraft-radiation --test trace
```

`output/furnace/paths.csv` has one row per point of a path (path,
emitter kind, point, what happens there, position, the fraction of the
energy left), for ParaView or any plotting tool.

[↑ back to top](#thermal-radiation)
