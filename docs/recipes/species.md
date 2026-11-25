# Species and reacting flows

Part of the [Recipes](../recipes.md): mixtures of ideal gases carried by the flow (steam furnace stage SF1). The keys are on the [key page](../cases/flow-keys.md#species) (`species.*`, `boundary.<name>.species`).

| # | Recipe |
|---|---|
| 278 | [Carry the species of a gas mixture with the flow](#species-transport) |

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
element C  inflow 2.130009e-7 kg/s, imbalance 1.385e-15 kg/s (6.50e-9 relative)
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
