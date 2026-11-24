# Porous regions

Part of the [Example gallery](../gallery.md): flow through porous matrices such as packed beds, foams and ceramic honeycombs, with a direction-dependent Darcy–Forchheimer resistance (D270, [porous plan](../assessments/porous-two-temperature-plan.md); two temperatures follow in PM2).

| # | Example | Stage | Feature |
|---|---|---|---|
| 273 | [A ceramic honeycomb in an air channel, and a porous block against the exact pressure drop](#porous-flow) | PM1 | `kind = "porous"`, `[regions.<name>.porous]`, `viscous`, `inertial`, `honeycomb`, `axes` |

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
