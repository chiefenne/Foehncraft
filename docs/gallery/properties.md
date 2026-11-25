# Expressions and materials

Part of the [Example gallery](../gallery.md): boundary values, body forces and viscosities given as expressions (D53), fluid properties from the material database and regions of cell zones (D102, D112), and walls with a mass flux (D124 W3).

| # | Example | Stage | Feature |
|---|---|---|---|
| 23 | [Expressions and materials: variable viscosity, blowing and suction, water from the database](#expressions-materials) | 5 | `[constants]`, `[tables]`, expressions, `[materials]`, `[regions]`, `foehn materials`, wall `mass_flux` |
| 277 | [Water and steam of IAPWS-IF97: the pressure-enthalpy chart, steam heated in a duct and superheated between hot plates](#steam-if97) | SF7a | `form = "if97"`, materials `steam` and `water-if97`, `foehn steam`, `materials.steam` |

---

<a id="expressions-materials"></a>
## 23. Expressions and materials: variable viscosity, blowing and suction, water from the database

<img src="../figures/gallery/expressions_materials.svg" alt="Left: plane Couette flow with the viscosity μ0 (1 + 4 y) from an expression. Middle: Couette flow with blowing through the lower wall and suction through the upper, Re = 4. Right: the viscosity of water in the material database, with the point of poiseuille_water.toml" width="800" height="298">

**A viscosity that varies in space.** `couette_viscosity.toml` gives the
viscosity as the string `"mu0 * (1 + a * y)"` with the constants
`mu0 = 0.1` and `a = 4`. The shear stress of plane Couette flow is the
same on every layer, so the velocity follows the viscosity:
u = U ln(1 + a y)/ln(1 + a), steeper near the lower wall where the fluid
is thinner (left; the dashed line is the constant-viscosity profile). The
viscosity is evaluated on every face; the observed order of the velocity
error is 1.83–1.96 on three meshes.

**Walls with a mass flux.** `couette_transpiration.toml` blows
0.2 kg/(m² s) through the lower wall and sucks the same through the upper
(`mass_flux = 0.2` and `-0.2`, D124 W3). The cross flow V = 0.2 m/s
pushes the shear layer against the upper wall:
u = U (e^(Re η) − 1)/(e^Re − 1), Re = V H/ν = 4 (middle). Observed order
1.86–1.96. The same key carries the evaporating mass of a wall film and the
species flux of a catalytic wall in later stages. A closed volume whose
walls add mass is refused before the run (BC026).

**Water from the database.** `poiseuille_water.toml` names
`material = "water"` at `temperature = 313.15` K: the density
992.2 kg/m³ and the viscosity 0.6527 mPa s come from the tables of the
built-in database (IAPWS-95 and the IAPWS 2008 viscosity formulation,
monotone cubic in T; right). `foehn materials` lists the database and
evaluates an entry at any temperature. The run writes the resolved
property program, with its sources, into `case.resolved.toml` and its
hash into every output file.

```
target/release/foehn flow examples/cases/flow/couette_viscosity.toml
target/release/foehn flow examples/cases/flow/couette_transpiration.toml
target/release/foehn flow examples/cases/flow/poiseuille_water.toml
target/release/foehn materials
target/release/foehn materials water --temperature 313.15
# the chart: the middle meshes with a line through the cell centres
for c in couette_viscosity couette_transpiration; do
  target/release/foehn flow examples/cases/flow/$c.toml --set 'verification.levels=[1]' \
      --set 'output.lines=[{name="mid",from=[1.0,0.03125,0.05],to=[1.0,0.96875,0.05],points=16}]'
done
python3 -m tools.plot.property_chart docs/figures/gallery/expressions_materials.svg
```

**Expressions against the built-in profiles.** Three more cases check
the expressions against the runs they replace: `poiseuille_expression.toml`
writes the inlet parabola as `"6 * u_mean * (y / h) * (1 - y / h)"` and
drives half of the flow by a body force expression; the velocity errors
equal those of `poiseuille_hex.toml` on all three meshes.
`taylor_green_expression.toml` prescribes the decaying vortex on the four
sides as expressions of x, y and t, evaluated at each time step (error
5.5e-4 at t = 0.4 s, as with `profile = "analytic"`).
`channel_ramp.toml` ramps the inflow with a table `ramp(t)`.

| Case | What the expression does | Check |
|---|---|---|
| `couette_viscosity.toml` | μ(y) | order 1.83–1.96 |
| `couette_transpiration.toml` | wall mass flux, blowing and suction | order 1.86–1.96 |
| `poiseuille_water.toml` | water at 313.15 K from the database, a region of the zone `fluid` | error 4.0e-4 |
| `poiseuille_expression.toml` | inlet profile and body force | errors of `poiseuille_hex.toml`, order 1.89–1.98 |
| `taylor_green_expression.toml` | velocity on the sides, in x, y and t | error 5.5e-4 |
| `channel_ramp.toml` | inflow ramped by a table in t | error 1.6e-3 |

How to use them: [recipe 31](../recipes/properties.md#expressions),
[recipe 32](../recipes/properties.md#materials).

---

<a id="steam-if97"></a>
## 277. Water and steam of IAPWS-IF97: the pressure-enthalpy chart, steam heated in a duct and superheated between hot plates

<img src="../figures/gallery/steam_ph_chart.png" alt="The pressure-enthalpy chart of water and steam from 1 kPa to 100 MPa and 50 to 4200 kJ/kg, coloured by the density on a logarithmic scale from 0.005 kg/m³ (blue, low-pressure steam) to 1000 kg/m³ (red, liquid): the saturation dome in black with the lines of quality 0.25, 0.5 and 0.75 dotted inside, its top at the critical point near 2090 kJ/kg and 22 MPa; white isotherms from 300 K to 1000 K run steeply through the liquid, flat across the dome and bend in the superheated steam; three white points at 4 MPa mark the inlets of the two steam cases and the outlet of the heated duct" width="800" height="517">

**The states.** `foehn steam` gives one state of water or steam from
(p, T), (p, h), (p, s), (p, x) or (T, x), or the saturation line at p or
T (tool `materials.steam`, MCP and Python through the registry). The
thermodynamic properties are the Industrial Formulation 1997 of IAPWS in
all five regions (Wagner et al. 2000; liquid, steam, the near-critical
region 3, the saturation line and steam up to 2273 K), the viscosity is
IAPWS R12-08 and the conductivity IAPWS R15-11 with its critical
enhancement, the surface tension IAPWS R1-76. Every verification point
of the releases is met to the rounding of their tables (5·10⁻⁹): regions
1, 2, 3 and 5 in v, h, u, s, c_p and w, the saturation line in both
directions, the boundary B23, eleven viscosities and the conductivities
away from the critical point. Region 3 is solved for the density from
its own Helmholtz equation, and (p, h) and (p, s) from the forward
equations, so a state is consistent with the basic equations to 10⁻¹²
(no backward equations). Two-phase states carry the equilibrium quality,
the homogeneous void fraction and density and three mixture viscosities.
The chart above is drawn from 16 500 such calls. The same body runs on
the GPU (entry `steam_state`), bitwise the CPU result.

<img src="../figures/gallery/steam_superheater_channel.png" alt="Left top: the temperature in a plane channel 60 mm long and 10 mm high, steam entering at 540 K on the left (blue) between two plates at 760 K: red layers grow from both walls and meet after about 20 mm, the core leaves at about 630 K (green), white streamlines straight along the channel. Left bottom: the density of the same run, 18.8 kg/m³ in the cold core (red) down to 11.8 kg/m³ at the walls (blue), the cold tongue shrinking along the channel. Right: the IF97 properties of steam at 4 MPa from the saturation temperature 523.5 K to 800 K, divided by their saturated-vapour values: viscosity rising to 1.7, conductivity dipping then rising to 1.45, density falling to 0.56, c_p falling steeply to 0.57 and the Prandtl number to 0.68" width="800" height="373">

**Steam heated in a duct (exact check).** `steam_heated_duct.toml` heats
superheated steam of the database material `steam` at 4 MPa in a duct
with slip walls from 600 K; the source adds 358.481 kJ/kg, which by IF97
is the step from h(4 MPa, 600 K) to h(4 MPa, 750 K). The energy equation
integrates the IF97 c_p, whose integral at constant pressure is the IF97
enthalpy: the steam leaves at 750.002 K (the 2 mK are the 0.25 K
sampling of the property curves), the mass flux is the same at both ends
to 10⁻¹², and the pressure the inlet needs, 20.55465 Pa, is the momentum
flux difference G²(1/ρ_out − 1/ρ_in) = 20.55466 Pa of the expanding
steam.

**A superheater in miniature.** `steam_superheater_channel.toml` lets
steam at 4 MPa and 540 K, 16.5 K above saturation, flow at 6 mm/s
between two plates at 760 K (Re 62, Pr 1.21). Density, c_p, μ and k
follow the temperature through IF97 with the variable density of stage
7c: the steam expands from 18.8 to 11.8 kg/m³ and speeds up, and leaves
at a mixed mean of 681 K. The run converges in 416 outer iterations and
the heat entering through the walls leaves with the steam to 4·10⁻⁸.
`ranks_water_if97_channel.toml` heats water of `water-if97` at 10 MPa
with IF97 c_p, μ and k on 2, 3 and 4 ranks: bitwise the one-rank
solution, in 121 iterations on every rank count.

```
cargo build --release -p foehncraft-cli
target/release/foehn steam --pressure 4e6 --temperature 700
target/release/foehn flow examples/cases/flow/steam_heated_duct.toml --out output/steam_duct
target/release/foehn flow examples/cases/flow/steam_superheater_channel.toml --out output/superheater
target/release/foehn flow examples/cases/flow/ranks_water_if97_channel.toml --out output/water_ranks
pip install matplotlib h5py numpy
python3 -m tools.plot.steam_states docs/figures/gallery output/superheater
```
