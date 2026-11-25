# Lagrangian sprays

Part of the [Example gallery](../gallery.md): droplets as Lagrangian parcels in the gas flow (L12, D240 to D246).

| # | Example | Stage | Feature |
|---|---|---|---|
| 242 | [Water droplets sprayed into an air duct](#spray-duct) | SP1 | `[spray]`, `[[injector]]`, `[liquid]`, drag of Schiller and Naumann, stick walls, mass book |
| 265 | [Where the liquid goes: datums and droplet statistics](#spray-datums) | SP1b | `[[datum]]` (cut plane, plane on faces, disc, STL surface, boundary, box), Rosin-Rammler sizes, `datums.csv`, histograms, D32 and Dv50, penetration, `series.json` |
| 280 | [A water spray that drives the air](#spray-two-way) | SP2 | `[spray.coupling] mode = "two-way"`, value-plus-slope momentum sources, passes inside the time step, coupling lag, momentum ledger pair, health numbers |

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

<a id="spray-two-way"></a>
## 280. A water spray that drives the air

<img src="../figures/gallery/spray_two_way.png" alt="Left: air speed in a room of 1.2 m by 1.5 m after 2 s, with white arrows of the air velocity and the droplets as dots coloured by their diameter; the falling spray in the middle has pulled the air into a jet of up to 6 m/s that spreads over the floor and rises along the side walls in two large eddies. Right, over time: the downward force of the spray on the air, which rises to about 0.8 N in the first half second and stays there; the parcel passes of each time step, mostly three; and the coupling lag of the last pass, mostly below the tolerance of 10⁻³." width="1200" height="660">

`spray_two_way.toml` is a room of still air, 1.2 m wide and 1.5 m high
(2D, 60 × 75 cells, SST k-ω, Δt = 4 ms). A nozzle under the ceiling sprays
50 g/s of water downwards at 15 m/s in a cone of 12° half angle, with
Rosin-Rammler sizes from 50 to 400 µm (X = 150 µm, n = 2.5). With
`[spray.coupling] mode = "two-way"` the drag on the droplets is a
momentum source of the gas (D247 to D252): the falling spray pulls the air
down into a jet of up to 6 m/s, which spreads over the floor and comes back
up along the side walls to the vents in the ceiling. One way (SP1) the air
would stay at rest.

**How the coupling works.** A parcel pass tracks every parcel through one
frozen gas state; each sub-step, which ends at every cell face, writes a
record of the momentum the gas receives and its slope with respect to the
gas velocity of the cell. The records are sorted by cell, parcel and
sub-step and summed in that order, so the sources do not depend on the
storage order of the parcels or the number of threads. The gas takes them
as value plus slope, S = E − β (u − u*), like a porous sink: the slope sits
on the diagonal and in the pressure correction, so a dense spray cannot
push the gas past the droplets, however large the step. The passes run
inside the outer iterations of the PISO step until the coupling lag is
below 10⁻³ (at most three, the default). After the last gas solve a
closing step gives each parcel the correction the gas received, so both
books hold the same numbers.

In this run the spray hands the air a steady 0.77 N downwards after the
first half second. Of the 100 g injected, 85.7 g land on the floor
(`wall_mode = "stick"`), 14.2 g are in the air and 0.1 g leave through the
vents; no parcel is lost. The steps take 2.96 passes on average; 46 of the
500 steps stop at the third pass with the lag above the tolerance (at most
3.1 × 10⁻³; SPR006 warns once and `spray.csv` counts them). The strongest
coupling of a cell, β Δt / (ρ V), is 0.77; the liquid fills at most 0.5 %
of a cell, with at most 16 parcels in it. The run takes about 40 s on
four cores.

**What is checked** (V28.1 to V28.4). A cloud shot into one cell of gas
at rest, with liquid loadings of 0.1 to 10 and time steps of 0.1 to 100
droplet response times: gas and droplets reach one common velocity, the
mixture keeps the momentum it was given, the force on the gas never
changes sign, and the momentum the gas received plus the momentum of the
parcels equals the momentum injected to round-off. One fixed pass is
bounded and conservative too. Against a reference of the two-mass
relaxation integrated by RK4 with small steps, the error falls at order
0.99 and 0.96 as Δt halves from 2 ms to 0.5 ms. The sources are bitwise
equal with the parcels stored in reverse order and on one thread (hex,
tetrahedral and polyhedral meshes), whole runs on 1 and 4 threads; the
slopes of one parcel add up to at most its mass. A parcel on a cell face
with the gas of both cells pushing it back across the face is not lost.
In the ledger the exchange is a pair (`spray drag` on the gas momentum,
`from the gas` on the parcel momentum) whose residual is zero to the bit,
and the parcel momentum account closes with the weight, the injection, the
walls and the parcels that left. The gas momentum accounts hold only this
term for now, so their imbalance is the force of the walls and the
pressure on the gas.

**What this version does not do.** Momentum only: no heat exchange (SP2b),
no turbulent dispersion (SP2c), no evaporation (SP5). Not with the energy
equation, conjugate heat transfer or an algorithm other than PISO (SPR003; the
steady spray mode is SP2s). The slope is the diagonal of the coupling: a
parcel that crosses cells in one step couples the neighbour cells only
through the passes. First order in time (the gas state at the new time
level, D248). One process. No time-step study tool yet (D254 (4)).

Files: `spray.csv` (passes, coupling lag, gas force, power and the health
numbers per step), `result.json` (coupling summary and ledger).
Picture: `python3 -m tools.plot.spray_two_way OUT_DIR output/spray_two_way`.
[Recipe 280](../recipes/sprays.md#spray-two-way),
[case row](../cases/flow.md#lagrangian-sprays).
