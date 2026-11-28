# Two fluids (volume of fluid)

Part of the [Example gallery](../gallery.md): two immiscible fluids such as water and air, separated by a sharp interface that moves with the flow (L3, volume of fluid). The volume fraction α of the first fluid is advected by THINC (default), CICSAM or upwind. Gravity and surface tension act as balanced face forces, and the pressure rows enforce the continuity of the volume. With the energy equation the liquid can evaporate and the vapour condense at the interface (MP2). Dispersed bubbles smaller than the cells take the two-fluid (Euler–Euler) model instead, with a momentum equation per phase (MP3).

| # | Example | Stage | Feature |
|---|---|---|---|
| 317 | [A water column collapses in air, against Martin & Moyce](#vof-dam-break) | MP1 | `[vof]`, `[[vof.phase]]`, `[[vof.region]]`, `scheme = "thinc"`, `surface_tension`, `vof.csv` |
| 358 | [Water evaporates and a vapour film boils, against exact solutions and Berenson](#vof-phase-change) | MP2 | `[energy]` with two phases, `[vof.phase_change]`, `phase_change_rate`, `examples/meshes/sphere_sector.py` |
| 368 | [Bubbly up-flow in a pipe with the two-fluid model, against Wang et al. (1987)](#two-fluid-pipe) | MP3 | `[two_fluid]`, `[two_fluid.gas]`, `[two_fluid.inflow.<boundary>]`, `volume_fraction`, `gas_velocity_x/y/z`, line columns `alpha_gas` |

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

---

<a id="vof-phase-change"></a>
## 358. Water evaporates and a vapour film boils, against exact solutions and Berenson

<img src="../figures/gallery/vof_phase_change_checks.png" alt="Three panels against time, each with the exact curve as a black line and the computed values as red circles on it, and the deviation below. Left, the Stefan problem: the vapour layer on a wall 10 K above saturation grows from 0.1 to 0.5 mm in 0.69 s, within 1 % of the exact layer after a dip of 3.4 % in the first steps. Middle, the sucking problem: the layer grows from 0.28 to 0.62 mm into water 2 K superheated, 1.4 to 1.8 % ahead. Right, Scriven's bubble: the radius grows from 0.1 to 0.3 mm in 8 ms in water 1.25 K superheated, 3.2 % ahead at 2 ms and 1.9 % at 8 ms." width="1400" height="680">

Three problems with exact solutions check the evaporation of water at 1 atm, where the vapour is 1604 times lighter than the liquid:
- **The Stefan problem** (`vof_stefan.toml`): a vapour layer between a wall 10 K above T_sat and saturated water grows by the heat it conducts, δ = 2λ√(α_v t). The vapour is at rest, the water leaves through the outlet.
- **The sucking problem** (`vof_sucking.toml`): the layer is at T_sat and the water 2 K superheated gives its heat to the interface, δ = 2β√(α_l t). The water carries its thermal layer, only 40 to 90 μm thick, towards the interface: the test of the energy convection across it.
- **Scriven's bubble** (`vof_scriven.toml`): a bubble grows in water 1.25 K superheated, R = 2β√(α_l t). The mesh is a 2° sector of the sphere with 1 μm cells, whose side planes pass through the centre.

How the phase change works: the interface cells are held at T_sat by an implicit sink, and the heat the sink takes, over the latent heat, is the rate. The vapour volume goes into the vapour cells beside the interface and the liquid sink into the liquid ones, so the continuity makes room for the vapour (D351).

What the runs give:
- **Stefan:** within 1.1 % of the exact layer after the start; the liquid volume closes with the outflow and the evaporated mass to 10⁻⁶.
- **Sucking:** 1.4 to 1.8 % ahead.
- **Scriven:** 3.2 % ahead at 2 ms, 1.9 % at 8 ms. The first steps, while the interface cells take their temperature, carry most of the lead.

The tests `phase_change_stefan`, `phase_change_sucking` and `phase_change_scriven` run shorter versions and hold 2 %, 2.5 % and 4 %.

<img src="../figures/gallery/vof_film_boiling.png" alt="Top: the speed around a boiling vapour film at four times, in a domain 79 mm wide and 236 mm high, with the interface as a white line. At 0.30 s the film on the wall is still thin with a wave that is thick at the two side walls. At 0.55 s two half bubbles rise from the film at the side walls with liquid flowing around them at up to 0.5 m/s. At 0.75 s and 1.00 s the first bubbles have pinched off and rise, a second pair follows, and a liquid jet runs up the middle between them. Bottom: the wall Nusselt number over the second of boiling: 3.9 at the start, a dip to 1.1 while the film thickens, then oscillating between 1.7 and 3.0 with each bubble release, mean 2.25 after 0.4 s, against the Berenson value of 2.62." width="1400" height="840">

`vof_film_boiling.toml` is film boiling on a horizontal wall 5 K above T_sat, in the set-up of Welch & Wilson (2000) with the near-critical fluid of the film boiling literature (density ratio 40, σ = 0.1 N/m, h_lv = 10 kJ/kg). The domain is one most dangerous Taylor wavelength wide, λ_d = 78.7 mm, and three high, with symmetry at the sides. The film starts with a cosine wave and the conduction profile across it.

What the run gives (80 × 240 cells, 1 s in 5000 steps of 0.2 ms, about 85 minutes on one core):
- **The instability:** the film thickens for 0.35 s, then the crests at the side walls rise as half bubbles, pinch off and leave; a new pair follows about every 0.22 s.
- **Heat transfer:** the wall Nusselt number Nu = q λ₀/(λ_v ΔT), λ₀ = √(σ/(g Δρ)), oscillates between 1.7 and 3.0 with the releases. Its mean after 0.4 s is 2.25, 14 % below the correlation of Berenson (1961), 2.62. One mesh only: no mesh study yet.

The colours are the speed, as the liquid and the vapour sit at T_sat outside the thin film.

```sh
target/release/foehn flow examples/cases/flow/vof_stefan.toml --out output/vof_stefan
target/release/foehn flow examples/cases/flow/vof_sucking.toml --out output/vof_sucking
target/release/foehn flow examples/cases/flow/vof_scriven.toml --out output/vof_scriven
target/release/foehn flow examples/cases/flow/vof_film_boiling.toml --out output/vof_film_boiling
python3 -m tools.plot.vof_phase_change checks docs/figures/gallery output/vof_stefan output/vof_sucking output/vof_scriven
python3 -m tools.plot.vof_phase_change film docs/figures/gallery output/vof_film_boiling
```

---

<a id="two-fluid-pipe"></a>
## 368. Bubbly up-flow in a pipe with the two-fluid model, against Wang et al. (1987)

<img src="../figures/gallery/two_fluid_pipe.png" alt="Left, two strips of the gas fraction in a vertical pipe 42 diameters long, the radius stretched and mirrored about the axis, flow from left to right. Top, j_G = 0.1 m/s: the gas enters uniform at 0.07, leaves the core within 7 diameters and gathers in a red band at the wall with up to 0.32; the core stays dark blue at about 0.02. Bottom, j_G = 0.4 m/s: the gas enters at 0.24, the core settles at about 0.11 after 25 diameters with a light wave on the way, and the wall band reaches 0.9. Right, the profiles at L/D = 35 from the axis to the wall against the measured points: at j_G = 0.4 the core lies on the measured 0.12 and the wall peak sits at the measured position but reaches 0.93 against 0.42; at j_G = 0.1 the core is 0.02 against 0.04 and the wall peak 0.32 against 0.16." width="1303" height="675">

Air bubbles of 3 mm rise with water in a vertical pipe of 57.15 mm (Wang, Lee, Jones & Lahey 1987): j_L = 1.08 m/s, j_G = 0.1 and 0.4 m/s, gas fraction measured across the radius at L/D = 35 (`bubbly_pipe_wang_jg01.toml`, `bubbly_pipe_wang_jg04.toml`). The two-fluid model solves one momentum equation for the water and one for the gas; Ishii–Zuber drag, Tomiyama lift, Hosokawa wall lubrication, Burns dispersion and Sato viscosity with k-ω SST for the water. The pipe is a 2° wedge, 240 × 24 cells with 1.2 mm across, both phases enter uniform with the turbulence of a developed pipe flow.

The lift drives the small bubbles towards the wall, where the wall force holds them about a bubble radius off: the wall peak of bubbly up-flow.

What the runs give (2000 steps of 2 ms, 4 s, about 11 minutes each on two cores):
- **j_G = 0.4 m/s:** core 0.114 against the measured 0.118 (−3 %); the wall peak at r/R = 0.965 where it was measured (0.953), but 0.93 against 0.42.
- **j_G = 0.1 m/s:** core 0.022 against 0.036 (−39 %); the wall peak at r/R = 0.965 (measured 0.96), 0.32 against 0.16.
- The steps converge in 7 to 9 outer iterations; the gas and water volume close with the boundaries.

The wall peaks are twice as high as measured: too much gas is held in the wall layer. The closures are the common set of one bubble size; the bubble-induced turbulence acts only as Sato's viscosity, with no source in k and ω, and the measured flows carry a spread of sizes. Size groups follow in MP5 to MP8 (D373 to D376). On a mesh graded to 0.3 mm wall cells, ten times smaller than a bubble, the gas piles into the wall cells: the averaged model wants cells of about half a bubble or more.

`two_fluid_column.toml` checks the drag exactly: in a uniform bubble column the slip is U_T √(1 − α) = 0.2247 m/s in every row (Ishii & Zuber, distorted regime), and each step converges in two outer iterations.

```sh
target/release/foehn flow examples/cases/flow/bubbly_pipe_wang_jg01.toml --out output/bubbly_pipe/jg01
target/release/foehn flow examples/cases/flow/bubbly_pipe_wang_jg04.toml --out output/bubbly_pipe/jg04
python3 examples/studies/bubbly_pipe/compare.py output/bubbly_pipe
python3 -m tools.plot.two_fluid pipe docs/figures/gallery output/bubbly_pipe
```

