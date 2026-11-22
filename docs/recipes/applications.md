# Engineering gallery

Part of the [Recipes](../recipes.md): run an application of the [Engineering Gallery](../applications.md) and draw its figures.

| # | Recipe |
|---|---|
| 161 | [Run the 2 L stirred tank of Fyfe et al.: power number and mixing time](#fyfe-stirred-tank-recipe) |

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
