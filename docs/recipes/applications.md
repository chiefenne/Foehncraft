# Engineering gallery

Part of the [Recipes](../recipes.md): run an application of the [Engineering Gallery](../applications.md) and draw its figures.

| # | Recipe |
|---|---|
| 161 | [Run the 2 L stirred tank of Fyfe et al.: power number and mixing time](#fyfe-stirred-tank-recipe) |
| 203 | [Mesh and run a multi-region heat exchanger from a STEP assembly](#shell-tube-hx-recipe) |

---

<a id="fyfe-stirred-tank-recipe"></a>
## 161. Run the 2 L stirred tank of Fyfe et al.: power number and mixing time

The case and its results are [application 161](../applications/rotating-machinery.md#fyfe-stirred-tank); the keys (`[[mrf]]`, `[turbulence]`, `[scalars]`, `[time]`, boundaries) are in the [flow case-file keys](../cases/flow-keys.md), the cases in the [case catalogue](../cases/applications.md).

```
pip install gmsh h5py scipy matplotlib
cargo build --release -p foehncraft-cli
python3 examples/applications/stirred_tank_fyfe/mesh.py --level 0
target/release/foehn flow examples/applications/stirred_tank_fyfe/case.toml --out output/fyfe
```

About 45 minutes on four cores for 2 500 iterations (194 827 cells). Look at `output/fyfe/boundaries.csv` for the impeller torque; N_P = 2π N |T| / (ρ N³ D⁵) with N = 200/60 1/s and D = 0.045 m, and compare with the band of 1.0 to 1.7 (measured values of a four-blade 45° pitched-blade turbine in a baffled tank, about 1.2 to 1.4). The torque oscillates by a few per cent (frozen rotor), so average the last iterations in `residuals.csv`.

The mixing time is a second, transient run that starts from the steady flow and carries a tracer released at the liquid surface:

```
target/release/foehn flow examples/applications/stirred_tank_fyfe/tracer.toml --out output/fyfe_tracer \
    --set 'run.restart="output/fyfe/restart.h5"'
python3 examples/applications/stirred_tank_fyfe/mixing_time.py output/fyfe_tracer
```

The script reads the `fields_<step>.vtkhdf` files (one every 2 s), prints θ₉₅ and writes `mixing_time.json` and `fyfe_mixing_time.png` with the measured value of the paper (18.7 s, Table 3, 200 rpm baffled) and the Grenville correlation (12.9 s at N_P 1.3). Look at the log after the first 100 steps: the run is expected to take one to two hours, so stop it and read the diagnosis rather than letting it run on. For another speed change `rpm` and `spin` in both cases; the reference times of the script are for 200 rpm.

The figure of the geometry is drawn by `python3 examples/applications/stirred_tank_fyfe/geometry_figure.py docs/figures/applications/fyfe_geometry.png`.

[↑ back to top](#engineering-gallery)

---

<a id="shell-tube-hx-recipe"></a>
## 203. Mesh and run a multi-region heat exchanger from a STEP assembly

The mesh and its checks are [application 203](../applications/heat-exchangers.md#shell-tube-hx), the case in the [case catalogue](../cases/applications.md).

```
pip install gmsh meshio pyvista
cargo build --release -p foehncraft-cli
python3 examples/applications/shell_tube_hx/mesh.py --level 2      # 983 249 cells, about 2 min
target/release/foehn mesh examples/applications/shell_tube_hx/mesh.toml
python3 examples/applications/shell_tube_hx/images.py examples/applications/shell_tube_hx/hx_L2.msh docs/figures/applications
```

`--level 0` writes a 0.48 M cell mesh for a quick test, `--level 3` a finer one (each level divides the cell size by 1.3). The script prints the volume of every zone and the area of every interface; compare the volumes with the CAD solids before meshing a new assembly.

**Your own assembly.** The steps carry over to any exchanger or cooled housing: close every opening with a plug just outside it (its flat face becomes the inlet or outlet), fragment a box around everything together with all solids, give each piece its zone from the solids that contain it, keep the box pieces inside the housing as fluid and drop the rest. Name every zone and every outer face as a physical group; leave the interfaces unnamed, the solver finds them from the zones. If a Boolean with a STEP solid gives pieces of the wrong volume, rebuild that part from OpenCASCADE primitives with the same dimensions, as `mesh.py` does for the tubes.

In the flow case each fluid zone becomes a `[regions.<name>]` of kind `fluid`, the metal zones one or more solid regions:

```toml
[mesh]
file = "hx_L2.msh"

[regions.metal]
kind = "solid"
zones = ["tubes", "tubesheets", "baffles", "housing"]
material = "steel"
```

The exchanger's own flow cases are `case_coflow.toml` and `case_counterflow.toml`: air on the shell side and water on the tube side, each fluid region with its own material (D194), the water turned round between them. The outlet pressures are absolute (`1.5e5` for the air, `1.0e5` for the water, D239), and `[regions.air.initial] temperature = "inlet"` starts the air at its 300 °C (D238). The case files name level 2; the results of application 203 are from level 0:

```
python3 examples/applications/shell_tube_hx/mesh.py --level 0      # 477 921 cells
target/release/foehn flow examples/applications/shell_tube_hx/case_coflow.toml \
    --set 'mesh.file="hx_L0.msh"' --threads 4 --out output/hx_coflow                  # about 90 min
target/release/foehn flow examples/applications/shell_tube_hx/case_counterflow.toml \
    --set 'mesh.file="hx_L0.msh"' --threads 4 --out output/hx_counterflow
pip install h5py pyvista pillow
xvfb-run -a python3 examples/applications/shell_tube_hx/flow_images.py \
    output/hx_coflow/fields.vtkhdf docs/figures/applications
```

`flow_images.py` writes the three flow images and the two images of the steel (`hx_steel_cutaway.png`, `hx_steel_bundle.png`, coloured from 20 to 30 °C; change `STEEL` in the script for another range). `flow.json` holds the heat flow and the mixed-mean temperature of every boundary (`solves[0].scalars`, the scalar `temperature`): the heat the air gives up is the sum of the enthalpy flows of the two shell nozzles, the water's the sum over the two tube nozzles. To continue a run, pass `--restart output/hx_coflow/restart.h5` with a larger `solver.max_iterations`; judge it by these heat flows and the air outlet temperature over the last few hundred iterations, not by the residuals, which stall near 10⁻³.

The renderings of the application page and the README come from ParaView (6.2, its own `pvpython`); the script reads the level 2 counter-flow result:

```
/Applications/ParaView-6.2.0.app/Contents/bin/pvpython \
    examples/applications/shell_tube_hx/render_images.py \
    output/hx_counterflow/fields.vtkhdf docs/figures/applications assets/readme
```

About 35 s on an iMac. It writes `hx_render_hero.jpg`, `hx_render_cut.jpg`, `hx_render_airpath.jpg` and `hx_render_baffle.jpg`, and with the third argument the README teaser `heat_exchanger_renderings.webp`, the hero and the air path cross-fading every 3.5 s (made by plain `python3`, whose Pillow needs WebP support). The steel is shaded as metal (physically based, metallic and roughness) under a studio environment map that the script draws itself; OSPRay path tracing is not used, because in ParaView 6.2 on macOS its path tracer loses the red channel. Each colour range is the 2nd to 98th percentile of what the picture shows. Change the cameras at the end of each block (`CameraPosition`, `CameraFocalPoint`) for other views and the seed count of `air_streamlines` for more or fewer streamlines.

