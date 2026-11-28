# Water and steam circuits (`examples/cases/circuits/`)

Part of the [Case catalogue](../cases.md): one-dimensional water and steam tube circuits of a boiler (stage SF7b, crate `foehncraft-circuits`). A circuit case is a TOML file of its own, run by the tool `circuits.run` (`foehn call circuits.run --input <file>.json`, inputs in [`examples/cases/tools/`](tools.md)). It writes `result.json` (nodes, tubes, drums, Newton history, warnings) and `profiles.csv` (every segment) into its output folder.

| Case | What it is | Discussed in |
|---|---|---|
| `parallel_tubes_laminar.toml` | cold water split between two unheated tubes of 2 m and 4 m in laminar flow: the Hagen–Poiseuille split 2:1 (exact) | [gallery 312](../gallery/circuits.md#water-steam-circuits) |
| `heated_riser_single_phase.toml` | water at 16 MPa heated by 200 kW/m² in a 10 m riser: the heat balance, the liquid column and the acceleration (exact) | [gallery 312](../gallery/circuits.md#water-steam-circuits) |
| `evaporator_tube.toml` | a once-through evaporator tube at 7 MPa, 400 kW/m², subcooled inlet to x = 0.47: subcooled and saturated boiling, Chen, Friedel, drift flux, Bowring CHF | [gallery 312](../gallery/circuits.md#water-steam-circuits), [recipe 312](../recipes/circuits.md#water-steam-circuits) |
| `wall_slab_tube.toml` | the water side of `flow/circuit_wall_slab.toml`: one tube at 10 MPa, 450 K inlet, 0.05 kg/s, heated by the slab above it (SF7c) | [gallery 356](../gallery/circuits.md#circuit-walls) |
| `wall_slab_boiling.toml` | the same tube at 7 MPa with 0.02 kg/s, 5 K subcooled: boiled by the slab to x = 0.43 (SF7c) | [gallery 356](../gallery/circuits.md#circuit-walls) |
| `evaporator_pass_drum.toml` | the drum loop at 7 MPa behind the two walls of `flow/circuit_evaporator_pass.toml`: two riser banks of 12 tubes per metre on one downcomer (SF7c) | [gallery 356](../gallery/circuits.md#circuit-walls), [recipe 356](../recipes/circuits.md#circuit-walls) |
| `drum_natural_circulation.toml` | a drum boiler loop at 7 MPa: downcomer, 30 heated risers, relief pipe; natural circulation ratio 16 | [gallery 312](../gallery/circuits.md#water-steam-circuits), [recipe 312](../recipes/circuits.md#water-steam-circuits) |

<a id="circuit-case-file-keys"></a>
## Circuit case-file keys

All values in SI units (Pa, K, J/kg, kg/s, m, W/m², W).

| Key | Meaning |
|---|---|
| `[circuit] name` | name of the case (default: the file name) |
| `[circuit] friction` | two-phase friction: `friedel` (default, Friedel 1979) or `homogeneous` (Churchill friction factor with the McAdams viscosity) |
| `[circuit] void` | `drift-flux` (default, Zuber–Findlay with Rouhani–Axelsson 1970) or `homogeneous` |
| `[circuit] boiling` | flow boiling heat transfer: `chen` (default, Chen 1966 with Butterworth's subcooled form) or `gungor-winterton` (1986) |
| `[circuit] subcooled` | vapour of subcooled boiling: `saha-zuber` (default, onset by Saha & Zuber 1974, true quality by Levy 1967) or `none` (thermal equilibrium) |
| `[circuit] chf` | critical heat flux: `bowring` (default, Bowring 1972), `table` (a look-up table in the layout of Groeneveld et al. 2007) or `none` |
| `[circuit] chf_table` | the CSV file of `chf = "table"`, relative to the case: header `p_kpa,g,x,chf_kw`, every grid point once, empty CHF cells allowed; the values are for an 8 mm tube and scaled by (0.008/D)^½ |
| `[circuit] chf_margin` | smallest CHF ratio q_CHF/q without the warning `CIR004` (default 1.3) |
| `[circuit] segments` | segments of each tube (default 40) |
| `[circuit] tolerance`, `max_iterations` | Newton tolerance on the scaled residuals (default 1e-10; pressures in 10⁵ Pa, enthalpies in 10⁵ J/kg, flows in the flow scale of the case) and iterations (default 100) |
| `[[nodes]] name`, `kind` | a node: `junction` (header: pressure and enthalpy from the balances), `pressure` (fixed pressure, the water it delivers has `temperature` or `enthalpy`) or `drum` (fixed pressure; the steam of the inflow leaves, the same mass of feedwater mixes with the separated water) |
| `[[nodes]] pressure` | pressure of a `pressure` or `drum` node; the start value of a junction |
| `[[nodes]] temperature`, `enthalpy` | the water a pressure node delivers, or the inflow `mass_flow` of a junction; a temperature is turned into an enthalpy at the node's pressure (a junction without one: the mean of the fixed pressures) |
| `[[nodes]] mass_flow` | mass flow into the circuit at a junction (negative: out of it) |
| `[[nodes]] feed_temperature` | feedwater temperature of a drum |
| `[[nodes]] elevation` | height of the node, for the report only (the tubes carry the rise) |
| `[[tubes]] name`, `from`, `to` | a tube; positive flow runs from `from` to `to` (a reversed flow is solved too) |
| `[[tubes]] length`, `diameter`, `rise` | length, inner diameter and the rise from `from` to `to` (−length ≤ rise ≤ length) |
| `[[tubes]] roughness` | absolute roughness (default 0) |
| `[[tubes]] count` | identical tubes in parallel (default 1); flows are per tube |
| `[[tubes]] segments` | segments of this tube (default `[circuit] segments`) |
| `[[tubes]] heat_flux` or `heat` | uniform heat flux on the inner wall, or the heat of one tube |
| `[[tubes]] heat_shape` | axial shape of the heat flux, pairs `[s, factor]` from s = 0 at `from` to 1 at `to`, linear between, scaled to a mean of one |
| `[[tubes]] unheated_inlet` | unheated length at the `from` end (default 0) |
| `[[tubes]] loss` | inlet loss coefficient K, Δp = K G²/(2ρ) |
| `[[tubes]] start_flow` | start value of the flow of one tube (default G = 1000 kg/(m² s)); give it near the expected flow for a natural circulation loop |
| `[output] profiles` | write `profiles.csv` (default true) |
