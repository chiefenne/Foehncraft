# Lagrangian sprays

Part of the [Example gallery](../gallery.md): droplets as Lagrangian parcels in the gas flow (L12, D240 to D246).

| # | Example | Stage | Feature |
|---|---|---|---|
| 242 | [Water droplets sprayed into an air duct](#spray-duct) | SP1 | `[spray]`, `[[injector]]`, `[liquid]`, drag of Schiller and Naumann, stick walls, mass book |
| 265 | [Where the liquid goes: datums and droplet statistics](#spray-datums) | SP1b | `[[datum]]` (cut plane, plane on faces, disc, STL surface, boundary, box), Rosin-Rammler sizes, `datums.csv`, histograms, D32 and Dv50, penetration, `series.json` |

---

<a id="spray-duct"></a>
## 242. Water droplets sprayed into an air duct

<img src="../figures/gallery/spray_duct.png" alt="Top: air speed in a duct of 2 m by 0.4 m with a nozzle in the upper wall at the left, and the droplets as dots coloured by their speed; the jet bends towards the outlet and the larger part lands on the lower wall. Bottom: the mass book of the run, injected, in the flow, stuck on the wall and left through the outlet, over time." width="800" height="465">

`spray_duct.toml` blows air at 4 m/s through a duct (SST k-ω, transient). A
nozzle in the upper wall sprays 10 g/s of water as droplets of 150 µm with a
speed of 10 m/s in a cone of 10° half angle. The droplets are parcels, each
standing for the same share of the mass flow; the injector sits in
`spray_duct_injectors.toml`, which the case includes (D246).

Each parcel follows the exact solution of its relaxation equation over every
sub-step: the drag uses the correlation of Schiller and Naumann (1935), the
sub-step ends at every cell crossing, and a small droplet that follows the
gas needs no small step (D245: the kernel is one source for CPU and GPU).
Droplets that reach the lower wall stay there (`wall_mode = "stick"`); the
others leave through the outlet. About 70 % of the liquid lands on the wall
within the 1.5 s of the case.

**What this version does not do.** The gas does not feel the droplets
(one way); there is no turbulent dispersion, no evaporation and no
breakup. Gas properties are sampled at the start of each sub-step, so the
tracking is first order in the gas state. One process only (D241). The
next stage (SP2) adds the coupling to the gas and the size distributions.

**What is checked.** The run books every kilogram: injected = in the flow +
stuck + left + lost, and stops with SPR005 when the book does not close
(`spray.csv` has it in every row). The unit tests compare the tracker with
a reference trajectory on hexahedral, tetrahedral, polyhedral and mixed
meshes, the Stokes terminal velocity and the exponential relaxation to the
gas velocity, and require identical bits on 1, 2 and 4 threads and after a
restart.

Files: `spray/parcels_<step>.vtkhdf` (ParaView), `spray.csv`; the parcels
are part of the restart file. [Recipe 242](../recipes/sprays.md#spray-duct),
[case row](../cases/flow.md#lagrangian-sprays).

---

<a id="spray-datums"></a>
## 265. Where the liquid goes: datums and droplet statistics

<img src="../figures/gallery/spray_datums.png" alt="Top: air speed in a duct of 2 m by 0.4 m with the droplets as dots coloured by their diameter; the large droplets fall in a curve onto the lower wall near the nozzle, the small ones fill the duct downstream; white lines mark a cut plane at x = 0.6 m, a disc at 0.9 m, a plane on the mesh faces at 1.2 m, a surface file at 1.6 m and a box under the nozzle. Bottom left: cumulative liquid through each datum over time, rising one after the other in the order of the flow. Bottom right: the mass distribution over the droplet diameter of the liquid that left through the outlet (small droplets, peak near 100 µm) and of the liquid on the floor (large droplets, 150 to 250 µm), with the Rosin-Rammler input of the nozzle." width="800" height="553">

`spray_datums.toml` is the duct of [242](#spray-duct) with a nozzle that
sprays 10 g/s of water with a Rosin-Rammler size distribution (X = 120 µm,
n = 2.5, from 10 to 400 µm; D242). Each parcel draws its diameter from the
mass distribution with the counter-based generator keyed by the parcel
identity (D121), so the sizes do not depend on threads or restarts.

Seven datums (`[[datum]]`, D262) ask where the liquid and the air go: a
plane cut through the cells at x = 0.6 m, the plane x = 1.2 m moved onto the
mesh faces, a disc of 0.1 m radius at x = 0.9 m, a surface file (STL) at
x = 1.6 m, the outlet, the lower wall and a box under the nozzle. Every
sub-step of a parcel that crosses a surface writes a crossing record; the
records are sorted by datum, parcel, sub-step and crossing and summed in
that order (D263), so the cumulative liquid through each datum has the same
bits on 1, 2 or 4 threads and after a restart.

After 1.5 s, of the 15 g injected, 8.3 g left through the outlet, 2.3 g
stick on the floor and 4.4 g are in the flow. The datum on the outlet is the
mass book's escaped liquid to the bit, the floor its deposited liquid. The
droplets sort by size: Dv50 is 91 µm at the outlet and 170 µm on the floor,
against 101 µm through the first plane (the input's Dv50 is 104 µm). The
Rosin-Rammler fit of the first plane's histogram gives X = 117 µm and
n = 2.6. In the flow the Sauter diameter D32 is 81 µm, and 95 % of the
liquid lies within 0.75 m of the nozzle along its axis (penetration, D265).

The gas side uses the same datums: through the cut plane the air mass flow
is 0.192 kg/s within 6 × 10⁻⁸ of the exact flow through the nearest face
layer (`cut_error`); the plane on the faces and the outlet are exact. The
disc (a 128-gon of equal area, of which the one-cell-thick 2D mesh holds the
middle strip) is cut through the cells, but its face layer is a staircase of
whole faces, so its `cut_error` is 4 %: here it measures the face layer, not
the cut. The STL surface is sampled finer than the cells and is within
0.1 % of the exact flow.

**What is checked** (V28.13 to V28.15). Planes, face layers, STL surfaces,
boundaries and boxes on hexahedral, tetrahedral and polyhedral meshes: the
liquid through every surface equals the liquid of the parcels that ended on
its far side, to round-off; the boundary datum equals the mass book. D32 of
the drawn sizes lies within four standard errors of the value of the
distribution, Dv50 within one bin; the sums are bitwise equal on 1, 2 and 4
threads and after a restart. In a flow run the outlet datum equals the
boundary mass flow bit for bit and the plane on the faces the inlet flow to
the solver tolerance.

**What this version does not do.** No measured anchor yet (V28.12, the
terminal velocity of Gunn and Kinzer, is still open); no enthalpy flow (the
energy equation with sprays comes with SP5); `datums.csv` only in transient
runs (steady runs report through `foehn call post.datums`); no size
classes, no time-averaged spray fields (D266, SP2) and no evaporation
numbers (SP5). Not with a sliding mesh. One process.

Files: `datums.csv`, `spray_hist_<datum>.csv`, `spray.csv` (penetration,
D32, centre of mass), `result.json` (all statistics with D20, D30, D43,
Dv10, Dv90, span, velocity statistics and the fit), `series.json`.
Picture: `python3 -m tools.plot.spray_datums OUT_DIR output/spray_datums`.
[Recipe 265](../recipes/sprays.md#spray-datums),
[case row](../cases/flow.md#lagrangian-sprays).
