# Engineering gallery

Part of the [Recipes](../recipes.md): run an application of the [Engineering Gallery](../applications.md) and draw its figures.

| # | Recipe |
|---|---|
| 161 | [Run the 2 L stirred tank of Fyfe et al. and compare it with the published CFD plane](#fyfe-stirred-tank-recipe) |

---

<a id="fyfe-stirred-tank-recipe"></a>
## 161. Run the 2 L stirred tank of Fyfe et al. and compare it with the published CFD plane

The case and its results are [application 161](../applications/rotating-machinery.md#fyfe-stirred-tank); the keys (`[[mrf]]`, `[turbulence]`, boundaries) are in the [flow case-file keys](../cases/flow-keys.md), the case in the [case catalogue](../cases/applications.md).

```
pip install gmsh h5py scipy matplotlib
cargo build --release -p foehncraft-cli
python3 examples/applications/stirred_tank_fyfe/mesh.py --level 0
target/release/foehn flow examples/applications/stirred_tank_fyfe/case.toml --out output/fyfe
```

About 45 minutes on four cores for 2 500 iterations (196 713 cells). Look at `output/fyfe/boundaries.csv` for the impeller torque; N_P = 2π N |T| / (ρ N³ D⁵) with N = 200/60 1/s and D = 0.045 m. The torque oscillates by about ±6 % (frozen rotor), so average the last iterations in `residuals.csv`.

To compare with the paper's CFD plane, download the plane files of the paper's supporting information (figshare doi:10.6084/m9.figshare.26405800, CC BY 4.0; not stored here) and run

```
python3 examples/applications/stirred_tank_fyfe/compare.py output/fyfe \
    "YZ plane 200rpm baffled Viscosity 0.0009.csv" --out output/fyfe
```

It writes `fyfe_speed_planes.png` and `fyfe_profiles.png`. The paper's plane has z = 0 at the bottom and z pointing down; the script turns it. For another speed, change `rpm` and `spin` in the case and take the plane file of that speed (the script's tip speed is for 200 rpm).
