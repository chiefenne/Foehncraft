- **Contact resistance, orthotropic solids and region heat sources
  (stage 7d2).** `[interfaces.<name>]` puts a contact resistance
  (`contact_resistance`, m²K/W) or a thin wall (`thin_wall` with
  `thickness` and `conductivity` or `material`, the resistance t/k) into
  the faces between two regions, fluid–solid or solid–solid; the face
  conductance is the series of the two half cells and the resistance, and
  the report gives the heat flow and temperature of every pair of regions
  (`interface-<solid>`, and `interface-<a>-<b>` between two solids). Solid
  materials may now have an orthotropic conductivity (three principal
  values along two axes, each a form of T), oriented per region with
  `[regions.<name>] axes`: the matrix takes the conductivity of the face
  normal, the tangential rest of the flux is a deferred correction.
  `[[energy.sources]]` heats one region with `power` (W, spread by volume)
  or `power_density` (W/m³). Codes CHT004 and ENE008.
- **Exact interface flux on skewed meshes.** The cells next to a region
  interface take one-sided least-squares gradients (the face temperature
  of their own side in place of the cell across), and the interface flux
  uses the face-normal temperatures of both sides (D157 §9.2): a
  piecewise-linear temperature is exact on any mesh and at any
  conductivity ratio. Check: the patch test `cht_orthotropic_patch.toml`
  (a laminate with axes turned by 30°, a fluid at rest and a block behind
  a contact resistance, on a mesh perturbed by 25 %) is exact to
  3·10⁻¹¹ K; V7.4 with a contact resistance and a thin wall matches the
  series formula to 10⁻⁹. Gallery and recipe 160: a 1 W chip on a thermal
  pad on a circuit board with copper planes (`cht_circuit_board.toml`).
- **Heat balance on skewed meshes.** The reported conduction through a
  boundary with a fixed temperature now includes the non-orthogonal part
  of the flux, as the matrix has it, so the balance of the energy equation
  closes on skewed meshes too.
- **Generated meshes.** With `perturb`, nodes on a face plane of a `zones`
  box now move only along the plane, so the interfaces between the zones
  stay planes.
