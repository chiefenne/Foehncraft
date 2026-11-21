# Rotating machinery

Part of the [Engineering Gallery](../applications.md): stirred tanks, pumps
and fans, with the rotor in a moving reference frame (D55) or on a sliding
mesh (D182). The features behind them are validated in the
[Example gallery](../gallery/flow.md) (Taylor–Couette flow, the baffled
mixer of gallery 100, the sliding-mesh studies).

| # | Application | Sector | Confidence |
|---|---|---|---|
| 131 | [Stirred tank with a three-blade hydrofoil: power and pumping number](#stirred-tank) | rotating machinery | demonstration |

---

<a id="stirred-tank"></a>
## 131. Stirred tank with a three-blade hydrofoil: power and pumping number

**Problem.** How much power a baffled stirred tank draws at a given speed,
and how much liquid the hydrofoil pumps down, which sets the mixing time and
the dead zones of the tank.

| | |
|---|---|
| Sector | rotating machinery |
| Geometry | T = 1 m, liquid height 1 m, dished base; impeller D = 0.4 m, three-blade hydrofoil; three baffles; sparger ring on three legs; feed pipe; source `examples/applications/stirred_tank/mesh.py` |
| Mesh | Gmsh tetrahedra, no prism layers; level 1: 607 032 cells, largest skewness 0.94 (12 faces above 0.9, D50 remediation on) |
| Fluid | water, ρ = 998.2 kg/m³, μ = 1.0016·10⁻³ Pa s; Re = ρND²/μ = 2.4·10⁵ |
| Boundary conditions | impeller and shaft: walls in the rotor frame (MRF, 90 rpm); vessel, baffles, sparger and feed pipe: walls; liquid surface: symmetry plane |
| Models | k-ω SST, automatic wall treatment; coupled solver |
| Run cost | 2 500 iterations, 2.8 h on one rank (2.1 GHz Intel Xeon) |
| Confidence | **demonstration**: C3 fails (torque still moves about ±4 %), C4 not run (one mesh) |
| Case files | `applications/stirred_tank/case.toml` (level 1); the smoke variant follows |

### Results

| Quantity | Value at 2 500 iterations | Mean over iterations 1 500–2 500 |
|---|---|---|
| Impeller torque | −1.05 N·m | −1.07 N·m |
| Power number N_P | 0.29 | 0.30 |
| Pumping number N_Q (downward, through the plane below the impeller) | 0.45 | – |
| Dead-volume fraction (speed below 5 % of the tip speed) | 13 % | – |
| Mean speed over tip speed | 0.13 | – |

### Images

Rendered from the level-1 run (800 iterations, kept in the project files;
the 2 500-iteration run is the one of the table, with the same flow
pattern); the script is `examples/applications/stirred_tank/images.py`.
The run is a frozen-rotor demonstration: the flow is not settled to the
check level, so read the images as the structure of the flow, not as
numbers.

<img src="../figures/applications/stirred_tank_speed_section.png" alt="Speed and streamlines in the vertical plane y = 0 of the stirred tank" width="800" height="826">

The impeller (centre) pumps the liquid down and outward along the dished
bottom; the liquid climbs along the wall and returns along the shaft, two
large recirculation loops on each side. The dark patch above the sparger is
the slow liquid between the sparger legs.

<img src="../figures/applications/stirred_tank_speed_horizontal.png" alt="Speed in the horizontal plane z = 0.22 m below the impeller of the stirred tank" width="800" height="723">

Below the impeller the discharge forms a ring of fast liquid around a slow
core; the three slow patches are the wakes between the baffles.

<img src="../figures/applications/stirred_tank_impeller.png" alt="Pressure on the three-blade hydrofoil impeller, with the hub, shaft, baffles and the sparger ring" width="800" height="738">

### Check card

| Check | Result |
|---|---|
| C1 Mass balance | pass: net imbalance 2·10⁻¹⁵ kg/s |
| C2 Energy balance | not run |
| C3 Iterative convergence | **fail**: torque moves 7 % over the last tenth of the iterations (rule: below 0.1 %); largest residual 2.8·10⁻³ |
| C4 Mesh sensitivity | not run (one mesh) |
| C5 Wall resolution | mean y⁺ 48 on the impeller, 72 on the vessel; the distribution is not yet checked |
| C6 Mesh quality | 12 faces above skewness 0.9, largest 0.94; remediation on |
| C7 Plausibility | the power number is in the range usually quoted for hydrofoil impellers; a plausibility check, not validation |
| C8 Independence of the domain | not run |
| C9 Health | MRF003: 5 676 cells where the rotor-interface flux does not close (the blades make an interface that is not a surface of revolution); RUN006: not converged. Both explained by the unsteady drift of the torque |

### What to trust, what not

- The mass balance and the power number are the most reliable numbers here.
  The power number carries the torque drift: about ±4 % from the mean.
- The pumping number comes from one plane below the impeller, so it is a
  flow estimate, not a field average.
- The dead-volume fraction depends on the mesh and on the threshold of 5 %
  of the tip speed. Treat it as indicative.

### Run it

From the repository root:

```
pip install gmsh h5py
python3 examples/applications/stirred_tank/mesh.py --level 1
target/release/foehn flow examples/applications/stirred_tank/case.toml
python3 examples/applications/stirred_tank/post.py output/case
```

### Limits and next steps

- The steady frozen rotor on this geometry does not settle; a longer
  averaging run, a sliding mesh, or a mesh-independence study (levels 2 and
  3) are the next steps.
- The images come from a 800-iteration run; they are re-rendered from the finished run by `images.py`.
