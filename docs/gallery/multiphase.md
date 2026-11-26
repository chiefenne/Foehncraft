# Two fluids (volume of fluid)

Part of the [Example gallery](../gallery.md): two immiscible fluids such as water and air, separated by a sharp interface that moves with the flow (L3, volume of fluid). The volume fraction α of the first fluid is advected by THINC (default), CICSAM or upwind. Gravity and surface tension act as balanced face forces, and the pressure rows enforce the continuity of the volume.

| # | Example | Stage | Feature |
|---|---|---|---|
| 317 | [A water column collapses in air, against Martin & Moyce](#vof-dam-break) | MP1 | `[vof]`, `[[vof.phase]]`, `[[vof.region]]`, `scheme = "thinc"`, `surface_tension`, `vof.csv` |

---

<a id="vof-dam-break"></a>
## 317. A water column collapses in air, against Martin & Moyce

<img src="../figures/gallery/vof_dam_break.png" alt="Left: four snapshots of the volume fraction of water in a tank 0.6 m long and 0.16 m high. At t = 0.05 s the column, 57 mm wide and 114 mm high at the left wall, starts to collapse. At t = 0.125 s it is a wedge reaching 0.17 m. At t = 0.225 s and t = 0.3 s it is a thin layer reaching 0.34 m and 0.48 m. Right: the surge front Z = z/a against T = t√(2g/a) on two meshes, a/16 and a/32, whose curves lie on top of each other, with the measured points of Martin & Moyce 1952 slightly below them." width="800" height="443">

`vof_dam_break.toml` is the classic test of Martin & Moyce (1952). A column of water a = 57.15 mm wide and 2a high stands at the left wall of a closed tank of air. At t = 0 it collapses under gravity. The surge front runs along the floor, and its position Z = z/a against T = t√(2g/a) is compared with the measured points of table 2 of the paper (n² = 2).

The run uses:
- **Fluids:** water and air, a density ratio of 828, with σ = 0.07 N/m.
- **Advection:** THINC with sub-steps at an outflow Courant number of at most 0.25.
- **Mesh and time:** 168 × 45 cells (h = a/16) and 1200 steps of 0.25 ms. It takes about 70 s on four cores.

What the run gives:
- **Agreement with the measurements:** the front follows them to within 5 % to 17 %, always a little ahead. Most of that is early, when the gate of the experiment was still lifting.
- **Mesh dependence:** halving the cells (a/32, 9 minutes) gives the same curve, so the gap is not a mesh effect.
- **Conservation:** the water volume is conserved to 2·10⁻¹⁴, and α stays in [0, 1] without clipping.

The test `dam_break_martin_moyce` repeats the run on h = a/8 and holds the front within 20 % of every measured point. `zalesak_disk_*` turns Zalesak's slotted disk once in a given rotation:

| Scheme | L1 error after one revolution |
|---|---|
| THINC (hex) | 0.048 |
| THINC (polygons) | 0.09 |
| CICSAM | 0.19 |
| Upwind | 1.10 |
| Semi-Lagrangian level set (computer-graphics benchmark) | 0.27 |

`vof_static_drop.toml` holds a drop of R = 2 mm at rest. It gives a Laplace jump of 32.2 Pa against the exact 35 Pa, with spurious currents of 0.18 m/s.

```sh
target/release/foehn flow examples/cases/flow/vof_dam_break.toml --out output/vof_dam_break
python3 -m tools.plot.vof_dam_break docs/figures/gallery output/vof_dam_break
```
