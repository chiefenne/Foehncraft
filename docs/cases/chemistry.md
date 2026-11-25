# Chemistry cases (`examples/cases/chemistry/`)

Part of the [Case catalogue](../cases.md): reactor case files of the gas-phase chemistry (steam furnace stage SF2). Run with `foehn chemistry reactor <file>` or `foehn chemistry equilibrium <file>`; each case checks its own result (`[reactor.check]`) and exits with 2 when a band is not met.

## Contents

- [Reactor cases](#reactor-cases)
- [Reactor case-file keys](#reactor-case-file-keys)

<a id="reactor-cases"></a>
## Reactor cases

| File | Mechanism | What it is | Discussed in |
|---|---|---|---|
| `chemistry/reactor_methane_air_gri30.toml` | GRI-Mech 3.0 (53 species, 325 reactions) | stoichiometric methane and air at 1 atm from 1400 K: ignition after 3.4 ms, end state 4·10⁻⁵ K and 4·10⁻⁹ (mole fraction) from the Gibbs equilibrium at 2697.7 K, element drift 2·10⁻¹³; 0.4 s | [gallery 289](../gallery/species.md#gas-chemistry), [recipe 289](../recipes/species.md#gas-chemistry) |
| `chemistry/reactor_methane_air_jl4.toml` | Jones & Lindstedt 1988 (4 steps) | the same mixture: ignition after 8 µs (a global mechanism has no induction chemistry), end state at the equilibrium of its eight species, 2816.9 K | [gallery 289](../gallery/species.md#gas-chemistry) |
| `chemistry/reactor_methane_air_wd2.toml` | Westbrook & Dryer 1981 (2 steps) | the same mixture: ends at 3164 K, 277 K above the equilibrium of its species, because its reverse CO step is an Arrhenius rate of its own and not the equilibrium constant (a flaw of the mechanism) | [gallery 289](../gallery/species.md#gas-chemistry) |

<a id="reactor-case-file-keys"></a>
## Reactor case-file keys

| Key | Type | Default | Meaning |
|---|---|---|---|
| `reactor.mechanism` | string | – | a built-in mechanism (`gri30`, `westbrook-dryer-2`, `jones-lindstedt-4`) or a YAML mechanism file (D103, the Cantera format), relative to the case file |
| `reactor.pressure` | Pa | 101325 | the constant pressure |
| `reactor.temperature` | K | – | the start temperature |
| `reactor.end_time` | s | – | the time to integrate |
| `reactor.mole_fractions` | table | – | the mixture by mole fractions (normalised; a species left out is 0) … |
| `reactor.mass_fractions` | table | – | … or by mass fractions … |
| `reactor.equivalence_ratio` | number | – | … or the fuel/oxidiser ratio over the stoichiometric one, with `fuel` and `oxidiser` (mole fractions) |
| `reactor.samples` | integer | 2000 | output times after the start |
| `reactor.spacing` | `linear`, `log` | `linear` | output times evenly in t, or evenly in log t from `end_time`/10⁶ |
| `reactor.integrator.rtol` | number | 1e-6 | relative tolerance of RODAS3 |
| `reactor.integrator.atol` | number | 1e-12 | absolute tolerance of the mass fractions |
| `reactor.integrator.atol_temperature` | K | 1e-6 | absolute tolerance of the temperature |
| `reactor.integrator.max_steps` | integer | 100000 | steps of one output interval before CHM002 |
| `reactor.check.temperature_difference` | K | – | the largest \|T_end − T_eq\| |
| `reactor.check.mole_fraction_difference` | number | – | the largest \|x_end − x_eq\| of any species |
| `reactor.check.equilibrium_temperature` | [K, K] | – | the range of the equilibrium temperature at the start enthalpy |
| `reactor.check.ignition_delay` | [s, s] | – | the range of the ignition delay (the time of the largest dT/dt) |
| `reactor.check.element_drift` | number | – | the largest relative change of an element mass |
