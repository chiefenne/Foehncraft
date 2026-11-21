# Energy equation and heat transfer

Part of the [Example gallery](../gallery.md): the energy equation of the fluid (stage 7, D150) with temperature-dependent c_p, k and μ, viscous heating and the thermal wall conditions, checked by manufactured and exact solutions ([specification](../requirements/stage7-energy-buoyancy-cht.md)).

| # | Example | Stage | Feature |
|---|---|---|---|
| 154 | [A heated cavity, conduction with k(T) and a Couette flow heated by its friction](#energy-fluid) | 7 | `[energy]`, `[boundary.<name>.energy]`, `[fluid] specific_heat`, `conductivity`, `viscous_heating`, `temperature_dependent`, `[solver.energy]` |
| 156 | [A heated wall cooled by a channel flow, and a composite wall](#cht-solids) | 7 | `[regions.<name>] kind = "solid"`, `[materials.<name>] kind = "solid"`, `[mesh.generate] zones`, `interface-<region>` |
| 158 | [A finned heat sink in a turbulent air channel, and a turbulent channel over a conducting wall](#cht-turbulent) | 7 | `[turbulence]` with solid regions, `[energy] wall_function`, `turbulent_prandtl` |

---

<a id="energy-fluid"></a>
## 154. A heated cavity, conduction with k(T) and a Couette flow heated by its friction

<img src="../figures/gallery/energy_cavity.png" alt="The temperature in the lid-driven cavity at Re = 100 with the lid at 310 K and the bottom at 300 K, Pr = 0.71, with the streamlines in white: the lid drags warm fluid down the right-hand wall and the primary vortex wraps it around its core" width="560" height="473">

**The heated cavity.** `cavity_heated.toml` is the Re = 100 cavity of
`cavity_re100.toml` with the lid held at 310 K, the bottom at 300 K and
adiabatic side walls; Pr = μ c_p/k = 0.71. The lid drags warm fluid down
the right-hand wall, the primary vortex wraps it around its core and the
heat leaves through the bottom. The run converges in 439 outer iterations
on 64 × 64 cells. The report gives the heat flow of each wall: 288.5 W/m²
enter at the lid and leave at the bottom, the mean Nusselt number is
Nu = q L/(k ΔT) = 2.05, and the heat balance closes to 7·10⁻¹² of the
largest flow. The temperature does not act on the flow yet (buoyancy is
PR 7b1).

<img src="../figures/gallery/energy_checks.png" alt="Left: the temperature of the manufactured conduction solution with k(T) = 1 + 0.01 (T − 300) on a perturbed 64 x 64 mesh. Middle: the conductivity that follows the temperature, from 1.0 at the walls to 1.5 in the centre. Right: the temperature across a Couette flow heated by its own friction on 8, 16 and 32 cells against the exact parabola" width="800" height="232">

**Conduction with k(T) (check V7.1).** `conduction_mms.toml` solves
∇·(k(T)∇T) + q''' = 0 in the unit square, perturbed by 15 %, with
k(T) = 1 + 0.01 (T − 300) W/(m K), the manufactured solution
T = 300 + 50 sin(πx) sin(πy) K and the source that belongs to it. The
conductivity is refreshed from the temperature in every outer iteration,
here from 1.0 at the walls to 1.5 in the centre. The L2 error is 0.112,
0.0258 and 0.00668 K on 16², 32² and 64² cells, observed orders 2.11 and
1.95 (the case requires 1.8 to 2.2); on the finest mesh 20.000 W of source leave through the four walls, with a
reported imbalance of 2·10⁻⁴ (the reported wall flux takes the orthogonal
part of the gradient, so the perturbed mesh shows a small gap).

**Couette flow heated by its friction (check V7.2).** In
`couette_viscous_heating.toml` the upper wall moves at U = 1 m/s over a
gap of h = 1 m (μ = 0.1 Pa s, k = 0.05 W/(m K)), both walls are held at
300 K and `viscous_heating = true` adds Φ = μ (∂u/∂y)². The temperature is
the parabola T = T₀ + μU²/(2k) η(1 − η), 0.25 K above the walls on the
centre line (Brinkman number 1):

| Cells across | L2 error of T (K) | Wall heat flux (W/m²) | Balance |
|---|---|---|---|
| 8 | 3.906e-3 | −0.0500 | 1e-9 |
| 16 | 9.766e-4 | −0.0500 | 3e-9 |
| 32 | 2.441e-4 | −0.0500 | 1e-8 |

The observed order is 2.000. Each wall takes half of the friction power
μU²/h = 0.1 W/m², and the conduction out of the walls equals the viscous
heating of the report.

**Properties that follow the temperature.** With a `material` the
specific heat, the conductivity and the viscosity come from the database
as functions of T (`temperature_dependent = true`, the default); μ(T)
reaches the momentum equation in every outer iteration. A Couette flow of
air between walls at 300 K and 600 K bends its velocity profile by 4 %
against the linear one, which comes back to below 10⁻⁶ with
`temperature_dependent = false` (`tests/energy.rs`). The convection
condition h (T∞ − T_w) gives the exact series resistance of a slab and a
film, and a `heat_flow` in W is spread over its boundary.

```
cargo build --release -p foehncraft-cli
target/release/foehn flow examples/cases/flow/cavity_heated.toml --out output/cavity_heated
target/release/foehn flow examples/cases/flow/conduction_mms.toml --out output/conduction \
    --set 'verification.levels=[2]' --set output.vtkhdf=true
for l in 0 1 2; do
  target/release/foehn flow examples/cases/flow/couette_viscous_heating.toml \
      --out output/couette_$l --set "verification.levels=[$l]" \
      --set 'output.lines=[{name="across",from=[1.0,0.0,0.05],to=[1.0,1.0,0.05],points=64}]'
done
python3 -m tools.plot.energy_fields docs/figures/gallery output/cavity_heated \
    output/conduction output/couette
```

---

<a id="cht-solids"></a>
## 156. A heated wall cooled by a channel flow, and a composite wall

<img src="../figures/gallery/cht_heated_wall.png" alt="Left: the temperature in a solid wall 0.2 m thick and the channel flow above it at Re = 100, with streamlines in white; the wall conducts ten times better than the fluid and is heated from below with 1000 W/m², the fluid warms from 300 K at the inlet towards the outlet. Right: the temperature along the top of the wall for conductivity ratios 0.1, 10 and 100; the better the wall conducts, the flatter and cooler its profile" width="800" height="243">

**A heated wall under a channel flow.** `cht_heated_wall.toml` is the
first conjugate case: one generated mesh of 80 × 30 cells holds a channel
1 m high and a solid wall 0.2 m thick below it, 4 m long, marked by a
`zones` box of the generator. The fluid enters at 300 K with 1 m/s
(Re = 100, Pr = 0.71); the outer face of the wall takes 1000 W/m². The
flow runs on the fluid cells, where the faces against the wall become the
no-slip boundary `interface-wall`; the temperature runs on all cells, with
the series conductance of the two half cells at every face between the
materials. All of the heater's 160 W (4 m × 0.04 m) cross the interface
and leave with the flow, and the heat balance of the whole domain closes
to 5·10⁻¹². With k_s/k_f = 10 (left) the heat spreads along the wall
towards the inlet, where the fluid is cold. The right-hand plot runs the
same case with three wall conductivities:

| k_s/k_f | Outer iterations | Heated face, mean (K) | Interface, mean (K) | Hottest cell (K) |
|---|---|---|---|---|
| 0.1 | 146 | 459.11 | 317.11 | 451.77 |
| 10 | 146 | 318.53 | 317.11 | 324.05 |
| 100 | 146 | 317.57 | 317.43 | 319.87 |

The interface mean is the same 317.1 K for the two lower ratios because
the fluid sees the same heat flow; a poor conductor (0.1) needs 142 K
across its thickness to pass it, a good one (100) spreads it so evenly
that the wall is almost isothermal.

**A composite wall (check V7.4).** `cht_composite_wall.toml` stacks a
solid with k = 200 W/(m K), a closed fluid layer at rest with k = 2 and a
solid with k = 20 000, 0.4 m each, between 300 K and 400 K: conductivity
ratios 1 : 100 : 10⁴. The temperature is piecewise linear and the series
heat flux is q = ΔT/Σ tᵢ/kᵢ = 495.0005 W/m². The run gives it at both
outer walls and at both interfaces to 10⁻⁹ relative, and the interface
temperatures 300.990 K and 399.990 K to 10⁻⁹ K (`tests/cht.rs`).

```
cargo build --release -p foehncraft-cli
target/release/foehn flow examples/cases/flow/cht_composite_wall.toml --out output/cht_wall
for r in 0.1 10 100; do
  k=$(python3 -c "print(14.084507042253522 * $r)")
  target/release/foehn flow examples/cases/flow/cht_heated_wall.toml --out output/cht_$r \
      --set "materials.wall-solid.conductivity=$k" --set output.vtkhdf=true
done
python3 -m tools.plot.cht_fields docs/figures/gallery output/cht
```

<a id="cht-turbulent"></a>
## 158. A finned heat sink in a turbulent air channel, and a turbulent channel over a conducting wall

<img src="../figures/gallery/cht_heat_sink.png" alt="Top: the temperature around an aluminium heat sink with three fins in a turbulent air channel, with streamlines in white; the metal is almost isothermal at about 329 K, warm air collects in the slow vortices between the fins and a long recirculation behind the last fin. Bottom: the speed, with the air accelerated to about 15 m/s over the fin tips and nearly at rest between the fins" width="800" height="440">

**A finned heat sink.** `cht_heat_sink_turbulent.toml` puts an aluminium
heat sink (k = 167 W/(m K)) on the floor of a plane channel 50 mm high:
a base 100 mm long and 5 mm thick with three fins 3 mm thick and 20 mm
tall. Air enters at 5 m/s and 300 K (k-ω SST with the automatic wall
treatment), and a chip under the base puts 2000 W/m² into the metal.
One mesh of 300 × 50 cells holds both: the flow runs on the air, the
temperature on everything, and every face between the air and the
aluminium takes the thermal wall function of Kader in series with the
half cell of the metal. The run converges in 1007 outer iterations
(2.5 min on 4 threads). All 2.000 W of the chip (per 10 mm of depth)
cross the interface and leave with the air, 0.68 K warmer on average;
the heat balance closes to 4·10⁻⁷. The aluminium conducts so well that
its top varies by only 0.7 K (328.5 to 329.2 K); the mean heat transfer
coefficient over the wetted surface is
q/(T_i − T_in) = 869.6 W/m² / 28.45 K = 30.6 W/(m² K). The fins block
half of the channel, so the air speeds up to 15 m/s over their tips,
while slow vortices fill the gaps between them, hold the warm air and
cool the inner fin faces poorly. The first-cell y⁺ on the metal runs
from 0.2 in the gaps to 24 at the fin tips, which the automatic wall
treatment and the thermal law bridge.

<img src="../figures/gallery/cht_turbulent_channel.png" alt="Left: the temperature of a turbulent channel flow over a solid wall held at 350 K on its outer face, the lowest 0.45 m of the 1.1 m high channel, stretched in height; the wall is warm and the fluid warms in a thin layer above it. Right: the cell temperatures across the wall and the fluid at x = 15 m for k-epsilon and k-omega SST: linear in the wall, a jump at the interface carried by the wall function, and the turbulent profile in the fluid" width="800" height="254">

**A turbulent channel over a conducting wall (check V7.36).**
`cht_turbulent_channel.toml` is half of a plane channel (H = 1 m, a
symmetry plane on top, Re = 2·10⁵, Pr = 0.71) over a solid wall 0.1 m
thick whose outer face is held at 350 K. The wall's resistance
t/k_s = 0.5 m²K/W is of the order of the fluid side's 1/h_w, so both
matter. The same channel without the solid, whose wall carries the
wall function in series with a convection condition h = k_s/t to 350 K,
is the independent path of the fluid-only code (stage 6f); the two must
agree, since the solid's axial conduction is negligible here
(`tests/cht.rs`):

| Model | Outer iterations | Heat flow, conjugate (W) | Heat flow, fluid only (W) | Interface, conjugate (K) | Wall, fluid only (K) |
|---|---|---|---|---|---|
| k-ε, scalable | 56 | 107.692 | 107.681 | 323.077 | 323.080 |
| k-ω SST, automatic | 64 | 106.023 | 106.006 | 323.494 | 323.498 |

The heat flow through the interface equals that of the outer face to
10⁻⁶, the drop across the wall is q t/k_s to 10⁻³, and the heat balance
closes to 10⁻¹⁰. On the right the temperature is linear in the wall and
jumps at the interface by the wall function's 1/h_w (first-cell y⁺ 70
to 85).

```
cargo build --release -p foehncraft-cli
target/release/foehn flow examples/cases/flow/cht_heat_sink_turbulent.toml --out output/cht_sink
target/release/foehn flow examples/cases/flow/cht_turbulent_channel.toml --out output/cht_channel_keps
target/release/foehn flow examples/cases/flow/cht_turbulent_channel.toml --out output/cht_channel_sst \
    --set turbulence.model=k-omega-sst --set turbulence.wall_treatment=automatic
python3 -m tools.plot.cht_turbulent docs/figures/gallery output/cht_sink \
    output/cht_channel_keps output/cht_channel_sst
```

