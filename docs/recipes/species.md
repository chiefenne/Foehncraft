# Species and reacting flows

Part of the [Recipes](../recipes.md): mixtures of ideal gases carried by the flow (steam furnace stage SF1). The keys are on the [key page](../cases/flow-keys.md#species) (`species.*`, `boundary.<name>.species`).

| # | Recipe |
|---|---|
| 278 | [Carry the species of a gas mixture with the flow](#species-transport) |
| 289 | [Run a mechanism in a constant-pressure reactor and find the equilibrium](#gas-chemistry) |

---

<a id="species-transport"></a>
## 278. Carry the species of a gas mixture with the flow

**Name the species.** The names are those of GRI-Mech 3.0 (53 species:
`H2`, `O2`, `H2O`, `CH4`, `CO`, `CO2`, `N2`, `AR`, `C2H6`, `C3H8`, `NO`,
`NO2`, …). The last one is the closure species, `1 − Σ` of the others;
take the one that dominates everywhere (N₂ in air).

```toml
[species]
names = ["CH4", "O2", "N2"]
initial = { mole_fractions = { O2 = 0.21, N2 = 0.79 } }
diffusion = "schmidt"     # rho D_i = mu/Sc_i; or "lewis": lambda/(c_p Le_i)
schmidt = 0.7             # or per species: { CH4 = 0.71, O2 = 0.76 }
turbulent_schmidt = 0.7   # with a turbulence model
properties = "mixture"    # c_p, lambda, mu of the mixture; "fluid": the fluid's
```

**Give every inlet its composition.** Mole fractions (numbers) or mass
fractions (numbers or expressions of x, y, z, t). A species left out is
zero; the closure species may be left out and takes the rest. An inlet
without a composition is an error (SPE001); walls and symmetry planes have
no species flux, outlets a zero gradient.

```toml
[boundary.fuel.species]
mole_fractions = { CH4 = 1.0 }
[boundary.air.species]
mass_fractions = { O2 = 0.233 }      # N2 takes 0.767
```

**Let the density follow the mixture.** With `[energy] density =
"variable"` the density is p₀ M/(R T) with the mixture's molar mass. The
fluid must be an ideal gas with its molar mass; its other properties are
only the start values when `properties = "mixture"`:

```toml
[materials.gas]
kind = "gas"
molar_mass = 0.0289647
density = { form = "ideal-gas" }
specific_heat = 1005.0
conductivity = 0.0263
viscosity = 1.85e-5
[fluid]
material = "gas"
[energy]
initial = 300.0
density = "variable"
```

Without `[energy]` set `properties = "fluid"`: an isothermal mixture with
the fluid's constant density and viscosity (a dilute tracer, or the
verification of the transport).

**Read the result.** The fields are `Y_<species>` for the transported
species (the closure species is 1 − Σ of them). The log gives the range of
every species and the closure of each element:

```
species CH4    Y in [-5.889788e-7, 1.000003e0]
element C  inflow 2.130009e-7 kg/s, imbalance 1.365e-15 kg/s (6.41e-9 of the inflow)
```

and `result.json` holds the accounts `species <name>` and `element <E>`
in `balance.last.ledger`, with the convected and the diffused mass per
boundary. An element that does not close to about the solver tolerance is
a sign of an unconverged run.

**Check the transport.** `[species] exact = { CH4 = "<expression>" }` gives
the error of a transported species as for a passive scalar, and
`[verification] scalar_max_error` makes it a check of the run:

    target/release/foehn flow examples/cases/flow/species_binary_diffusion.toml
    target/release/foehn flow examples/cases/flow/species_coaxial_jet.toml
    python3 -m tools.plot.species_fields output/img output/species_coaxial_jet output/species_binary_diffusion

Mixing only, no reactions yet (SF2, SF3); every species has a Schmidt or a
Lewis number (no mixture-averaged diffusion yet, SF2h); one fluid region.
[Gallery 278](../gallery/species.md#species-transport).

---

<a id="gas-chemistry"></a>
## 289. Run a mechanism in a constant-pressure reactor and find the equilibrium

**Write the reactor.** A reactor case file has one table `[reactor]`: the
mechanism, the pressure, the start temperature, the time and the mixture
(the keys are on the [chemistry case page](../cases/chemistry.md#reactor-case-file-keys)):

```toml
[reactor]
mechanism = "gri30"          # or "jones-lindstedt-4", "westbrook-dryer-2", "my_mech.yaml"
pressure = 101325.0          # Pa
temperature = 1400.0         # K
end_time = 0.1               # s
equivalence_ratio = 1.0      # or mole_fractions = { CH4 = 0.095, O2 = 0.19, N2 = 0.715 }
fuel = { CH4 = 1.0 }
oxidiser = { O2 = 0.21, N2 = 0.79 }
samples = 600
spacing = "log"              # output times evenly in log t
```

**Run it.** The trace (time, temperature, the mole fraction of every
species) goes to `output/<case>/reactor.csv`; the log gives the ignition
delay (the time of the largest dT/dt), the end state against the
equilibrium at the start enthalpy, the element drift and the steps:

    target/release/foehn chemistry reactor examples/cases/chemistry/reactor_methane_air_gri30.toml

```
reactor: gri30 (53 species, 325 reactions), p = 101325 Pa, T₀ = 1400 K, 0.1 s
  ignition delay 3.427905e-3 s
  T end 2697.701 K, equilibrium 2697.701 K, difference 4.105e-5 K
  largest |x_end − x_eq| 4.338e-9, element drift 2.229e-13
  2811 steps, 0 rejected, 0.379 s
```

**Make it a check.** `[reactor.check]` gives bands (`temperature_difference`,
`mole_fraction_difference`, `equilibrium_temperature`, `ignition_delay`,
`element_drift`); a band not met is reported and the exit code is 2, as for
a flow case with `[verification]`.

**Find the equilibrium.** The same file, at its enthalpy (the adiabatic
flame state) or with `--constant-temperature` at its temperature:

    target/release/foehn chemistry equilibrium examples/cases/chemistry/reactor_methane_air_gri30.toml

From an agent or a script the two are the tools `chemistry.reactor` and
`chemistry.equilibrium` ([tools.md](../tools.md)), with the `[reactor]`
table inline or a case file:

    target/release/foehn call chemistry.equilibrium --input examples/cases/tools/equilibrium_methane_air.json

**Bring your own mechanism.** A YAML mechanism file of the Cantera format
(D103) is read directly: `units`, the species of the first phase, NASA-7
thermo and Lennard-Jones transport records, elementary, three-body and
fall-off (Lindemann, Troe) reactions with efficiencies and `orders`. A
species the file does not define is taken from the built-in GRI-Mech 3.0
data, so a global mechanism needs only its reactions. PLOG, Chebyshev and
SRI rates and specific colliders are refused (CHM001), as is a reaction
that does not conserve its elements. If the integrator stops (CHM002),
raise `[reactor.integrator] max_steps` or loosen `rtol`.

**Mind the global mechanisms.** They are calibrated for flame speed and
heat release. Jones & Lindstedt ignites at once and Westbrook & Dryer
ends 277 K above its own equilibrium (its reverse CO step is not the
equilibrium constant): do not use them for ignition or for burnt-gas
temperatures. The images of the gallery:

    python3 -m tools.plot.reactor_traces output/img

[Gallery 289](../gallery/species.md#gas-chemistry).

