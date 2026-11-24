# Porous regions

Part of the [Example gallery](../gallery.md): flow through porous matrices such as packed beds, foams and ceramic honeycombs, with a direction-dependent Darcy–Forchheimer resistance (D270) and gas and matrix with temperatures of their own (D271, [porous plan](../assessments/porous-two-temperature-plan.md)).
Part of the [Example gallery](../gallery.md): flow through porous matrices such as packed beds, foams and ceramic honeycombs, with a direction-dependent Darcy–Forchheimer resistance (D270, [porous plan](../assessments/porous-two-temperature-plan.md); two temperatures follow in PM2).

| # | Example | Stage | Feature |
|---|---|---|---|
| 273 | [A ceramic honeycomb in an air channel, and a porous block against the exact pressure drop](#porous-flow) | PM1 | `kind = "porous"`, `[regions.<name>.porous]`, `viscous`, `inertial`, `honeycomb`, `axes` |
| 275 | [Charging a honeycomb heat store, and Schumann's bed against the exact solution](#porous-two-temperatures) | PM2 | `porous.solid`, `heat_transfer_coefficient`, `specific_area`, `solid_conductivity`, `honeycomb.nusselt`, field `matrix_temperature` |

---

<a id="porous-flow"></a>
## 273. A ceramic honeycomb in an air channel, and a porous block against the exact pressure drop

<img src="../figures/gallery/porous_honeycomb_channel.png" alt="Top: air speed with streamlines in a channel 60 cm long and 20 cm high with a honeycomb block of 15 cm by 12 cm in its middle, channels along the flow; the air inside the block runs straight and leaves only through its back face. Middle: the same block with the same resistance in every direction; the air that enters its front face turns out through the upper and lower faces. Bottom: the pressure on the axis of both runs, linear through the honeycomb, curved through the isotropic block." width="800" height="800">

`porous_honeycomb_channel.toml` blows air at 10 cm/s through a channel
with a block of the "Solarstone" ceramic honeycomb brick (square channels
of 6.9 mm on a pitch of 13.6 mm, open area ε = 0.256). The case gives the
matrix in one line, `honeycomb = { channel = 6.9e-3, pitch = 13.6e-3,
length = 0.15 }`: the resistance along the channels is that of laminar
square-duct flow, f_D Re = 56.9, the entrance loss of a brick 0.15 m long
comes from K(∞) = 1.43 (both Shah & London 1978), and across the channels
the resistance is 10⁴ times larger (`blocking`). The axes of the region put
the channels along x.

The block resists strongly, so most of the air goes around it in both
runs. The difference is inside: in the honeycomb the air that enters runs
straight through and leaves at the back face, and the pressure falls
linearly; with the same resistance in every direction (the second command,
`blocking = 1`) the air turns out of the block through its upper and lower
faces, and the pressure on the axis falls fast near the front face. The
velocity is the superficial one (volume flow over the whole cross section),
continuous across the block faces.

<img src="../figures/gallery/porous_block_convergence.png" alt="Left: the pressure along a duct with a porous block between x = 0.4 m and 0.6 m on 50, 100 and 200 cells against the exact piecewise linear profile. Right: the relative error of the pressure drop against the cell size, falling with the second power." width="800" height="293">

`porous_darcy_block.toml` checks the drop itself: a duct with slip walls
and a block of K⁻¹ = 10⁶ 1/m² and F = 10 1/m, through which
Δp = (μ K⁻¹ U + ½ ρ F U²) L = 4.8 Pa exactly. The resistance jumps at the
block faces, which the pressure gradient of the neighbouring cells sees
over one cell; the error of the drop is 0.41 % on 100 cells and falls with
the square of the cell size. A duct filled with the matrix has no jump and
gives the drop to the solver tolerance (10⁻⁸), along and across the
channels of a direction-dependent matrix.

**What this version does not do.** Axes turned against the mesh axes are
refused (POR002): the Rhie–Chow factor of the solver is one number per
cell and cannot hold the tensor, and a turned honeycomb stalls; a tensor
factor is follow-up work. Channels across the main flow (a block that
closes the channel) converge slowly. No scalars or turbulence in porous cases
(POR002; energy came with [275](#porous-two-temperatures)), one process, CPU only. The kinetic-energy ledger does not book the
closes the channel) converge slowly. No energy, scalars or turbulence in
porous cases yet (POR002; the two temperatures of gas and solid are stage
PM2), one process, CPU only. The kinetic-energy ledger does not book the
work of the porous sink.

**What is checked.** The exact drop of a filled duct, along and across the
channels of a direction-dependent matrix (to 10⁻⁸), the order of the drop
over a block on 50, 100 and 200 cells (observed order above 1.5, error
below 0.5 % on 100 cells), the honeycomb values against the formulas, a
turned tensor stays symmetric, and the checks of the keys (POR001, POR002).

Files: `fields.vtkhdf` (ParaView), `result.json`. [Recipe 273](../recipes/porous.md#porous-flow),
[case rows](../cases/flow.md#porous-regions).

---

<a id="porous-two-temperatures"></a>
## 275. Charging a honeycomb heat store, and Schumann's bed against the exact solution

<img src="../figures/gallery/porous_storage_charge.png" alt="Left four panels: a section of a heat store 40 cm wide and 1 m high with a stack of ceramic honeycomb bricks from 15 to 85 cm. The matrix temperature after 1500, 3000 and 6000 s: a horizontal hot band moves down from the top of the stack and reaches the middle after 6000 s. Fourth panel: the gas temperature after 3000 s with white streamlines, a swirl in the upper plenum, straight downward lines in the stack. Right: gas and matrix temperature on the axis at the three times; the matrix lags the gas by up to about 30 K at the front." width="800" height="406">

`porous_storage_charge.toml` charges a section of a honeycomb store: air
at 573 K enters through an opening in the top right corner at 0.4 m/s,
crosses a plenum, flows down through a stack of the "Solarstone" bricks
(6.9 mm channels on a pitch of 13.6 mm, ε = 0.256, channels vertical) and
leaves through an opening in the bottom left. Gas and bricks start at
293 K. The stack is a porous region with `solid = "ceramic"`, so every cell
of it carries two temperatures: the gas, with ε ρ c and ε k, and the
matrix, with (1 − ε) ρ_s c_s and the conductivity of the brick walls,
(1 − ε) k_s along the channels and (1 − √ε) k_s across (D272). The two
exchange heat with h a, h = Nu k / d_h with Nu = 2.98 (laminar square
ducts, Shah & London 1978) and a = 4 ε / d_h: h = 15 W/(m² K) and a = 148 m²/m³ here.

The honeycomb holds the air in its channels, so the flow in the stack is
straight down and uniform within 1 % whatever the plenum does, and the
front is flat. The bricks hold about 7 000 times more heat than the air
in the channels, so the thermocline moves 7 000 times slower than the air:
half way down the stack after 6 000 s (3.4 mm/min).
Ahead of the front the matrix lags the gas by up to about 30 K; behind
it the two are within a few kelvin.

<img src="../figures/gallery/porous_schumann.png" alt="Gas and matrix temperature along a bed 1 m long after 1000, 2000 and 3000 s, every third of 100 cells as circles and squares, on the solid and dashed lines of the exact solution of Schumann (1929): the gas falls from 400 K at the inlet, the matrix lags below it." width="800" height="404">

`porous_schumann.toml` checks the coupling against the exact solution of
Schumann (1929): air at 400 K into a bed at 300 K, 0.1 m/s, h a =
603 W/(m³ K), no conduction. Gas and matrix lie on the exact curves at
every time; the L1 error falls with the square of the cell size (8.2·10⁻⁴,
2.1·10⁻⁴, 5.2·10⁻⁵ of the step on 25, 50 and 100 cells). The heat balance
of every step closes to 10⁻¹³ of the largest flow; the matrix is an account
of its own (`bed matrix` in the storage).

**How it is built.** The matrix temperature is one more row per porous cell
in the one energy matrix of the run (D271): a monolithic system, not a
loop between two equations, so a large h a costs no extra iterations and
the two temperatures relax to each other exactly (V3). The storage of the
matrix is in enthalpy form with the time scheme of the run (Euler or BDF2).

**What this version does not do.** The matrix is adiabatic at the faces of
the stack: it exchanges heat with the gas only, not with a housing or a
neighbouring solid (D274, PM3). Gas properties are constant (no density
that follows T in porous cases, POR002), no radiation in the channels
(R1b), no cycle-averaged or accelerated energy modes, no frozen-flow march
yet (PM3). The extended system is solved on the host, one process. In the
store the swirl in the upper plenum does not converge within the 10 outer
iterations of a step at a Courant number of 2 400; the stack flow is not
affected. A CHT restart does not carry the matrix temperature yet (a flow
restart does).

**What is checked.** V3: a closed box of gas and matrix relaxes to the
exact discrete decay, heat kept to 10⁻⁹. V4: Schumann's bed on 25, 50 and
100 cells, error below 10⁻⁴ on 100 cells and observed order above 1.8. V5:
with a large h a, gas and matrix are within 0.1 K and the front moves at
the speed of the mixture within 0.01 m. Every step of these runs closes its
heat balance to 10⁻⁹.

Files: `fields.vtkhdf` (`temperature`, `matrix_temperature`), `line_axis.csv`,
`heat.csv`, `result.json`. [Recipe 275](../recipes/porous.md#porous-two-temperatures),
[case rows](../cases/flow.md#porous-regions).
