# Water and steam circuits

Part of the [Recipes](../recipes.md): the one-dimensional water and steam tube circuits of a boiler (stage SF7b, crate `foehncraft-circuits`). A circuit case is its own TOML file ([keys](../cases/circuits.md#circuit-case-file-keys)), run through the tool registry.

| # | Recipe |
|---|---|
| 312 | [Solve a boiler circuit: an evaporator tube and a drum loop with natural circulation](#water-steam-circuits) |
| 356 | [Put a circuit behind the walls of a flow case](#circuit-walls) |

---

<a id="water-steam-circuits"></a>
## 312. Solve a boiler circuit: an evaporator tube and a drum loop with natural circulation

1. Write the nodes and the tubes. A drum fixes the pressure level and
   separates the steam; junctions are headers; `count` puts identical
   risers in parallel:

   ```toml
   [circuit]
   boiling = "chen"           # or "gungor-winterton"
   chf = "bowring"            # or "table" with chf_table = "lut.csv"

   [[nodes]]
   name = "drum"
   kind = "drum"
   pressure = 7.0e6
   feed_temperature = 490.0

   [[nodes]]
   name = "lower_header"
   kind = "junction"

   [[tubes]]
   name = "downcomer"
   from = "drum"
   to = "lower_header"
   length = 24.0
   rise = -22.0
   diameter = 0.25
   loss = 1.5
   start_flow = 40.0          # near the expected circulation

   [[tubes]]
   name = "risers"
   from = "lower_header"
   to = "drum"
   count = 30
   length = 20.0
   rise = 20.0
   diameter = 0.05
   heat_flux = 1.5e5
   start_flow = 1.3
   ```

2. Run it with a JSON input that names the case and the output folder:

   ```sh
   cargo build --release -p foehncraft-cli
   target/release/foehn call circuits.run --input examples/cases/tools/circuits_drum_loop.json
   ```

   The drum's steam flow and circulation ratio, and every tube's flow,
   pressure drop and smallest CHF ratio go to standard error. The result
   envelope goes to standard output.

3. Read the results. `output/drum_natural_circulation/result.json` has the
   nodes (pressure, enthalpy, temperature; a drum with its steam flow,
   riser flow and circulation ratio) and the tubes (flow per tube and
   in total, mass flux, the friction, gravity and acceleration parts of
   the pressure drop, the outlet state, the hottest wall and the smallest
   CHF ratio with its place). `profiles.csv` has every segment: pressure,
   enthalpy, temperature, equilibrium and flowing quality, void, heat
   flux, wall temperature, heat transfer coefficient, regime, CHF and CHF
   ratio.

4. Draw the profiles:

   ```sh
   target/release/foehn call circuits.run --input examples/cases/tools/circuits_evaporator.json
   python3 -m tools.plot.circuits docs/figures/gallery output/evaporator_tube output/drum_natural_circulation
   ```

5. React to the warnings. `CIR004` means a segment's CHF ratio is below
   `chf_margin` (below 1 the wall is past the CHF: lower the heat flux,
   raise the flow or the pressure). `CIR002` means the Newton solve did not
   converge: give `start_flow` near the expected flows. `CIR003` means a
   tube ran out of the property range: too little flow for its heat.

---

<a id="circuit-walls"></a>
## 356. Put a circuit behind the walls of a flow case

1. Write the circuit as in recipe 312. The tubes that sit behind the
   walls get their heat from the flow. A `heat_flux` given for them is
   only the start of the coupling. Give start flows near the expected
   ones: at small heat a drum loop can also settle with the flow reversed
   (risers down, steam up the downcomer).
2. In a steady flow case with `[energy]`, name the circuit and map each
   wall boundary to one tube. `start` is the point of the tube's `from`
   end and `end` that of its `to` end. Every face of the boundary must
   project onto the axis between them:

   ```toml
   [circuits]
   file = "../circuits/evaporator_pass_drum.toml"
   every = 20          # outer iterations between circuit solves
   tolerance = 1e-4    # coupling residual Σ|Q_3D − Q_circuit|/Σ|Q_3D|

   [[circuits.walls]]
   boundary = "xmin"
   tube = "risers_left"
   start = [0.0, 0.0, 0.5]
   end = [0.0, 6.0, 0.5]
   wall_resistance = 1.5e-4   # metal between the face and the water, per face area
   fouling = 5.0e-5           # scale on the inner wall, per inner area
   ```

   Leave out `[boundary.xmin.energy]`, because the circuit gives the
   condition. A bank of `count` tubes shares the faces of its height. Set
   `wall_resistance = 0` when the tube metal is a solid region of the
   mesh.
3. Run it:

   ```sh
   foehn flow examples/cases/flow/circuit_evaporator_pass.toml --out output/evaporator_pass
   ```

   The log shows the frame "Water and steam circuit (SF7c)" with the
   circuit solves, the coupling residual, the drum steam and the coupled
   tubes. The run converges only when both the flow residuals and the
   coupling residual are below their tolerances; otherwise `CIR006`.
4. Read the outputs. `circuit.csv` holds one row per circuit solve.
   `circuit/result.json` and `circuit/profiles.csv` are the circuit as in
   recipe 312. `circuit/walls.csv` holds every mapped face with its
   segment, U, the water temperature, the heat flux into the water and the
   gas-side wall temperature. `result.json` carries the summary under
   `circuits`, and the ledger has the account of region `circuit` with
   the exchange pair 2901.
5. Draw the gallery image:

   ```sh
   pip install matplotlib h5py numpy
   python3 -m tools.plot.water_walls docs/figures/gallery output/evaporator_pass
   ```

Steady runs on one rank only, CPU or GPU. A transient case is refused
(`CIR005`; transient circuits are stage SF18), and so is a run on ranks.
