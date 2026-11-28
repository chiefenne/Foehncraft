# Water and steam circuits

Part of the [Example gallery](../gallery.md): the one-dimensional water and steam tube circuits of a boiler in the crate `foehncraft-circuits` (stage SF7b of the [steam furnace path](../assessments/steam-furnace-development-path.md)), on the IAPWS-IF97 states of SF7a. A circuit is a graph of tubes between headers, pressure nodes and steam drums. Its flows, pressures and enthalpies are solved together, including the natural circulation of a drum loop.

| # | Example | Stage | Feature |
|---|---|---|---|
| 312 | [Flow boiling in an evaporator tube, and the natural circulation of a drum boiler](#water-steam-circuits) | SF7b | `circuits.run`, `[[nodes]] kind = "drum"`, `[circuit] boiling`, `friction`, `void`, `subcooled`, `chf` |
| 356 | [A convective evaporator pass on a drum loop: the circuit behind the 3D walls](#circuit-walls) | SF7c | `[circuits]`, `[[circuits.walls]]`, `wall_resistance`, `fouling` |

---

<a id="water-steam-circuits"></a>
## 312. Flow boiling in an evaporator tube, and the natural circulation of a drum boiler

<img src="../figures/gallery/circuits_evaporator.png" alt="Left: a vertical tube 12 m high drawn as a band coloured by the void fraction, blue (no vapour) in the lowest 1.5 m, through green and yellow near 3 to 4 m to dark red (void 0.9) at the outlet. Right top: along the tube the equilibrium quality rises in a straight line from -0.1 to 0.47; the flowing quality of the Levy profile leaves zero at 1.5 m, before the equilibrium quality reaches zero at 2.1 m, and then joins it; the void fraction rises steeply to 0.6 at 4 m and 0.88 at the outlet. The background is green for subcooled boiling up to 2.1 m and sand for saturated boiling beyond. Right middle: the bulk temperature rises from 531 K to the saturation temperature 559 K at 2.1 m and stays there; the wall temperature lies 3 to 14 K above saturation, highest where the bulk reaches saturation, then falls slowly as convective boiling grows. Right bottom: the CHF ratio of Bowring on a log scale falls from 12 at the inlet to 3.6 at the outlet, above the dashed margin 1.3" width="800" height="485">

**The tube.** `evaporator_tube.toml` is a once-through evaporator tube,
12 m high with a 20 mm bore. Water enters at 530 K (29 K subcooled) at
7 MPa, 0.35 kg/s (G = 1114 kg/(m² s)), and the wall is heated by
400 kW/m². The tube is marched in 80 segments with the IF97 states (left
and upper right). Wall heat brings the energy to the equilibrium quality
x_eq = 0.47 at the outlet, exact to round-off. Vapour appears at the
Saha–Zuber point (Saha & Zuber 1974), 0.6 m before the bulk is
saturated. The flowing quality follows the profile of Levy (1967). The
void fraction comes from the drift flux of Zuber & Findlay with the C₀
and V_gj of Rouhani & Axelsson (1970).

**The wall.** The wall temperature (middle) comes from Gnielinski in the
liquid. Once the wall passes saturation it comes from Chen (1966) with
Butterworth's subcooled form, or from Gungor & Winterton (1986) with
`boiling = "gungor-winterton"`. Both correlations give a superheat of 3
to 14 K, within a factor of two of each other along the whole tube (a
test). The critical heat flux of Bowring (1972) stays 3.6 times above
the heat flux at the outlet (bottom). At 1 MW/m² and 0.25 kg/s the same
tube reaches it. The run then warns with `CIR004`, and the wall beyond
the CHF takes the vapour convection of Dougall & Rohsenow (1963), more
than 50 K hotter.

<img src="../figures/gallery/circuits_drum_loop.png" alt="Left: a drum boiler loop drawn to scale: a grey drum at 22 m, a downcomer 0.25 m wide down to the lower header at 0 m in blue (no vapour), five of the thirty risers up to the upper header at 20 m coloured by the void fraction from blue at the bottom through cyan near 10 m to green (0.5) at the top, and the relief pipe in green back to the drum; title: steam 7.29 kg/s, circulation ratio 16.1. Right top: along a riser the equilibrium quality rises from -0.015 to 0.06, the flowing quality leaves zero near 4 m and the void fraction rises from 4 m to 0.5 at 20 m. Right bottom: the bulk temperature rises from 555 K to the saturation temperature 560 K at 5.4 m; the wall temperature jumps from 555 K to 564 K where the heating starts at 1 m, peaks at 568 K at 5.5 m and falls to 565 K at the top" width="800" height="529">

**Natural circulation.** In `drum_natural_circulation.toml` a drum at
7 MPa feeds a downcomer 22 m down to a lower header. Thirty risers, each
20 m high with a 50 mm bore and heated by 150 kW/m² above their first
metre, carry wet steam to an upper header, and a relief pipe returns it
to the drum. Nothing sets the flow. It follows from the weight of the
cold downcomer column against the lighter two-phase column of the risers,
less the friction and the inlet losses. The drum takes the steam out of
the riser flow and mixes the same mass of feedwater at 490 K into the
water of the downcomer. The Newton solve of the whole loop converges in 5
iterations. It gives 117 kg/s around the loop and 7.29 kg/s of steam, a
circulation ratio of 16 (left). The mass, energy and pressure balances
close exactly (tests): the heat of the risers equals the steam flow
times (h″ − h_feed) to 10⁻⁶, and the pressure drops of the three tubes
add to zero to 10⁻³ Pa.

**Exact checks.** `parallel_tubes_laminar.toml` splits a laminar flow
between two tubes 2 m and 4 m long. They take 2/3 and 1/3 of the flow as
Hagen–Poiseuille demands, with the inlet pressure drop to 10⁻⁶.
`heated_riser_single_phase.toml` heats water at 16 MPa in a riser. The
outlet enthalpy is the heat balance to round-off, and the pressure drop
is the sum of friction, the liquid column and the acceleration
G²(1/ρ_out − 1/ρ_in).

```sh
cargo build --release -p foehncraft-cli
target/release/foehn call circuits.run --input examples/cases/tools/circuits_evaporator.json
target/release/foehn call circuits.run --input examples/cases/tools/circuits_drum_loop.json
python3 -m tools.plot.circuits docs/figures/gallery output/evaporator_tube output/drum_natural_circulation
```

**Gaps.** The CHF look-up table of Groeneveld et al. (2007) is not
shipped. `chf = "table"` reads it from a CSV file the user supplies, and
Bowring is the built-in default. Void and wall temperatures are not yet
compared with measurements (Bartolomei et al. 1982). Post-dryout heat
transfer is the plain Dougall–Rohsenow vapour convection, with no
non-equilibrium. Supercritical water uses Gnielinski with bulk
properties, with no correction for the pseudo-critical region. The
circuit is steady, with no flow instabilities or transients (SF18). The
coupling to the 3D furnace walls is SF7c.

---

<a id="circuit-walls"></a>
## 356. A convective evaporator pass on a drum loop: the circuit behind the 3D walls

<img src="../figures/gallery/circuit_evaporator_pass.png" alt="Left: a vertical gas pass 0.2 m wide and 6 m high, drawn eight times wider than its scale, coloured by the gas temperature: yellow-white (1600 K) at the bottom left where the hotter gas enters, orange (1400 to 1450 K) in the core all the way up, dark purple thermal boundary layers (1000 to 1100 K) along both walls. Beside each wall a bar for its bank of risers, coloured by the void fraction from blue (zero) at the bottom through cyan, green and yellow to red (0.09) at the top. Second panel: the heat flux into the water along the height, peaking at 36 kW/m² on the left wall and 23 kW/m² on the right wall near the inlet and falling to 16 and 14 kW/m² at the top. Third panel: the gas-side wall temperature from 568 K (left) and 565 K (right) at the bottom to 563 K at the top, the inner walls at 560 to 561.5 K, and the water at 559 to 559.4 K. Fourth panel: the void fraction rises from zero near 0.3 m to 0.095 (left) and 0.086 (right) at the top; ten times the equilibrium quality rises from -0.02 at the bottom through zero near 1.3 m. Title: steam 0.142 kg/s, circulation ratio 144; left bank 125 kW and 11.2 kg/s, right bank 95 kW and 9.2 kg/s; 8 circuit solves, coupling residual 2.1e-07" width="1430" height="836">

**The case.** `circuit_evaporator_pass.toml` is a convective pass of a
boiler: flue gas rises at 20 m/s between two water walls, 0.2 m apart and
6 m high (2D, 1 m deep, k-ε with the thermal wall function of Kader).
The gas enters hotter on the left (1600 K) than on the right (1200 K), as
behind a burner that fires off the middle. Each wall is the outer face of
a bank of 12 risers per metre (40 mm bore) of one drum loop at 7 MPa
(`../circuits/evaporator_pass_drum.toml`): a drum, a downcomer, a lower
and an upper header, the two banks and a relief pipe. Between the gas and
the water lie 6 mm of steel (`wall_resistance = 1.5e-4` m² K/W) and a thin
scale (`fouling = 5e-5` m² K/W).

**The coupling.** Every wall face is mapped to the segment of its tube by
its height. The water side enters the 3D energy equation as a Robin
condition q″ = U (T_water − T_b) on each face, inside the energy matrix
(D349). U is the series of the wall, the scale and the inner heat transfer
coefficient of the segment (Gnielinski, Chen), with the ratio of the gas
side to the inner area of the tubes. Every 20 outer iterations the heat of
the faces is summed per segment and handed to the circuit. Its Newton
solve starts from the last solution and takes one or two iterations. The
run converges in 143 outer iterations with 8 circuit solves, and the
coupling residual Σ|Q_3D − Q_circuit|/Σ|Q_3D| is 2·10⁻⁷. Explicit
exchange of a heat flux or a wall temperature, lagged by a solve, would
oscillate or need heavy relaxation (the flaw the [development
path](../assessments/steam-furnace-development-path.md) section 8 names).

**The result.** Neither bank's flow is prescribed. The two banks share
one downcomer, and each takes the flow that its heat drives: the hotter
left bank gets 125 kW and 11.2 kg/s, the right bank 95 kW and 9.2 kg/s.
The drum makes 0.142 kg/s of steam at a circulation ratio of 144. The
heat flux into the water peaks just above the inlet, where the thermal
boundary layers are thin (36 kW/m² left, 23 kW/m² right), and falls to
about 15 kW/m² at the top. Without the radiation of the gas (SF4) these
heat fluxes are those of convection alone, a tenth of what a furnace wall
sees, so the circulation ratio is high. The ledger books the circuit's
own thermal account: heat from the walls, feedwater in and steam out. It
closes to 5·10⁻¹² of the heat, and its exchange pair with the conduction
of the two walls leaves 0.05 W of 220 kW, the coupling residual.

**Exact check.** `circuit_wall_slab.toml` is a slab 2 m long and 50 mm
thick, held at 700 K on top, on a tube of water at 10 MPa that enters at
450 K. The coupled run converges in 12 iterations and 7 circuit solves.
The outlet is 519.3815 K. An independent march of the same heat balance
(RK4 on 4000 steps with the IF97 states and Gnielinski,
`tests/water_walls.rs`) gives 519.3816 K. Every face sends the heat flux
of the series resistance of the slab and the water side to 3·10⁻⁴; only
the axial conduction near the adiabatic ends of the slab departs from it.
On a boiling tube at 7 MPa (`../circuits/wall_slab_boiling.toml`) the
same slab boils water to x = 0.43 in 12 iterations. Most of the tube is
then at the saturation temperature, and the wall is that of Chen.
