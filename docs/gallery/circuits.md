# Water and steam circuits

Part of the [Example gallery](../gallery.md): the one-dimensional water and steam tube circuits of a boiler in the crate `foehncraft-circuits` (stage SF7b of the [steam furnace path](../assessments/steam-furnace-development-path.md)), on the IAPWS-IF97 states of SF7a. A circuit is a graph of tubes between headers, pressure nodes and steam drums. Its flows, pressures and enthalpies are solved together, including the natural circulation of a drum loop.

| # | Example | Stage | Feature |
|---|---|---|---|
| 312 | [Flow boiling in an evaporator tube, and the natural circulation of a drum boiler](#water-steam-circuits) | SF7b | `circuits.run`, `[[nodes]] kind = "drum"`, `[circuit] boiling`, `friction`, `void`, `subcooled`, `chf` |

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
