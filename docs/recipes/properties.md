# Expressions and materials

Part of the [Recipes](../recipes.md): boundary values, body forces and viscosities as expressions, and fluid properties from the material database. The keys are on the [key page](../cases/flow-keys.md#expressions-materials-and-regions-d53-d102-d112).

| # | Recipe |
|---|---|
| 31 | [Write boundary values, body forces and the viscosity as expressions](#expressions) |
| 32 | [Take the fluid properties from the material database](#materials) |
| 277 | [Use water and steam of IAPWS-IF97: states, saturation lines and the `if97` material form](#steam-if97) |

---

<a id="expressions"></a>
## 31. Write boundary values, body forces and the viscosity as expressions

**A value as a string.** Any velocity, pressure, mass flux, mass flow,
body force or viscosity may be a string instead of a number (D53). Name
the numbers once in `[constants]`:

```toml
[constants]
u_max = 1.5
h = 0.41
omega = "2 * pi * 4.0"        # constants may use other constants

[boundary.inlet]
type = "velocity-inlet"
velocity = ["4 * u_max * y * (h - y) / h^2", 0.0, 0.0]

[boundary.outlet]
type = "static-pressure-outlet"
pressure = "100 * sin(omega * t)"   # transient cases only

[flow]
body_force = [0.0, "-1.2 * 9.81", 0.0]

[fluid]
viscosity = "1e-3 * (1 + 4 * y)"    # x, y, z; not t
```

**Tables.** A `[tables.<name>]` is a function of one argument, piecewise
linear and constant beyond its ends; inline or from a two-column file:

```toml
[tables.ramp]
x = [0.0, 0.5]   # t (s)
y = [0.0, 1.0]

[tables.lift]
file = "lift_curve.csv"   # relative to the case file

[boundary.inlet]
type = "velocity-inlet"
velocity = ["6 * y * (1 - y) * ramp(t)", 0.0, 0.0]
```

**Check before you run.** `foehn flow --check` parses and compiles every string
and reports the mistakes with their place and the nearest name:

```
target/release/foehn flow --check examples/cases/flow/channel_ramp.toml
target/release/foehn flow examples/cases/flow/channel_ramp.toml
target/release/foehn flow examples/cases/flow/taylor_green_expression.toml
target/release/foehn flow examples/cases/flow/couette_viscosity.toml
```

| Code | Meaning |
|---|---|
| EXPR001 | syntax error (with a caret under the place) |
| EXPR002 | unknown or reserved name, cycle between constants (with "did you mean") |
| EXPR003 | a variable the key does not take: `t` in a steady case or in the viscosity, x, y, z in a mass flow |
| EXPR004 | a table that is missing or not increasing |
| EXPR005 | a value that is not finite where it is evaluated |

**Walls with a mass flux.** `mass_flux` on a wall (kg/(m² s), into the
domain) blows or sucks through it; the wall's tangential velocity stays
as given. The flow needs a pressure boundary that takes up the net mass
(else BC026):

```toml
[boundary.porous_wall]
type = "wall"
mass_flux = "0.2 * step(x, 0.5)"   # blowing downstream of x = 0.5 m
```

The reserved names: the variables `x y z t`, the constants `pi e`, the
functions, `if clamp and or not`, and the names of solution fields that
later stages bring into expressions (`u v w p T rho mu k omega epsilon
nut Y`).

<a id="materials"></a>
## 32. Take the fluid properties from the material database

**List and evaluate.** The database has dry air (ideal gas, Sutherland),
water (IAPWS tables), engine oil, ethylene glycol and glycerol, each with
the source of its data and its temperature range:

```
target/release/foehn materials
target/release/foehn materials water --temperature 313.15
target/release/foehn materials air --temperature 350 --pressure 2e5
target/release/foehn --json materials        # for scripts
```

**Use one.** Name it in `[fluid]`; until the energy equation (L1) the
properties are evaluated once, at `temperature` and `operating_pressure`:

```toml
[fluid]
material = "water"
temperature = 313.15          # K, default 293.15
operating_pressure = 101325   # Pa, the default
```

**Change a property.** A `[materials.<name>]` table with `base` takes a
database entry and overrides what it names; a new material needs no base:

```toml
[materials.warm-oil]
base = "engine-oil"
viscosity = { form = "power-law", value_ref = 0.25, t_ref = 320.0, exponent = -6.0, source = "own measurement" }

[materials.test-fluid]
kind = "liquid"
density = 1050.0
viscosity = { form = "table", temperature = [280.0, 300.0, 320.0], values = [4e-3, 2.5e-3, 1.7e-3], interpolation = "monotone-cubic" }
```

**Regions.** `[regions.<name>]` gives cell zones their material (D112):
fluid regions only until the energy equation (REG001); every zone in one
region (REG002). A region with another viscosity gives a viscosity per
face:

```toml
[regions.core]
kind = "fluid"
zones = ["fluid"]
material = "water"
```

```
target/release/foehn flow examples/cases/flow/poiseuille_water.toml
grep -A12 "resolved property program" output/poiseuille_water/case.resolved.toml
```

`case.resolved.toml` lists the resolved properties with their form,
temperature and source; `properties_hash` and `materials_database` in the
provenance of every output file identify the property program and the
database it came from.

---

<a id="steam-if97"></a>
## 277. Use water and steam of IAPWS-IF97: states, saturation lines and the `if97` material form

**One state.** Give the pressure with one other quantity, in SI units
(Pa absolute, K, J/kg, J/(kg K)):

```
foehn steam --pressure 4e6 --temperature 700        # superheated steam
foehn steam --pressure 4e6 --enthalpy 3.2e6         # T from (p, h)
foehn steam --pressure 7e6 --quality 0.25           # wet steam
foehn steam --temperature 373.15 --quality 0        # saturated liquid
foehn steam --pressure 101325                       # the saturation line
foehn --json steam --pressure 4e6 --entropy 6.5e3   # JSON for programs
```

The state lists the region, phase, ρ, h, u, s, c_p, c_v, the speed of
sound, the equilibrium quality (below 0 for a subcooled liquid, above 1
for superheated steam), the void fraction, μ, k, Pr and the surface
tension; a two-phase state adds the saturated phases and the mixture
viscosities of McAdams, Cicchitti and Dukler. From Python or an agent
the same call is the tool `materials.steam`:

```python
import foehncraft as fc
s = fc.call("materials.steam", pressure=4e6, temperature=700.0)["result"]["state"]
```

A state outside IF97 (above 100 MPa, below 273.15 K, above 1073.15 K
beyond 50 MPa or above 2273.15 K) is refused with MAT010.

**Steam or water as the fluid of a case.** The database has `steam`
(the vapour branch) and `water-if97` (the liquid branch); both evaluate
IF97 at the operating pressure of the case, which with absolute pressure
values comes from the outlet:

```toml
case_format = "0.7"

[fluid]
material = "steam"
temperature = 600.0          # the case temperature (K)

[energy]
density = "variable"         # ρ(T) from IF97; c_p, μ, k follow T anyway

[boundary.outlet]
type = "static-pressure-outlet"
pressure = 4.0e6             # absolute: IF97 is evaluated at 4 MPa
```

Any material can take the form for any of its four properties, and
`phase` picks the branch below the critical pressure (default: the
phase of the case temperature):

```toml
[materials.feedwater]
kind = "liquid"
density = { form = "if97", phase = "liquid" }
specific_heat = { form = "if97", phase = "liquid" }
viscosity = { form = "if97", phase = "liquid" }
conductivity = { form = "if97", phase = "liquid" }
```

The liquid branch ends at T_s(p) and the vapour branch starts there; a
temperature beyond the end of its branch takes the saturated state and
warns with MAT003. Boiling itself (a fluid crossing T_s) is not a
single-phase property: the 1D tube circuits of SF7b and the two-fluid
model of SF10 bring it. The cases are
[`steam_heated_duct.toml`](../../examples/cases/flow/steam_heated_duct.toml),
[`steam_superheater_channel.toml`](../../examples/cases/flow/steam_superheater_channel.toml)
and [`ranks_water_if97_channel.toml`](../../examples/cases/flow/ranks_water_if97_channel.toml)
([gallery 277](../gallery/properties.md#steam-if97)).
