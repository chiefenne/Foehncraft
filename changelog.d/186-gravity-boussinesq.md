- **Gravity and the Boussinesq buoyancy (stage 7b1).** `[gravity]` with
  `vector` and `reference_point`, and `[buoyancy]` with
  `model = "boussinesq"`, `expansion_coefficient` (default from the
  fluid's density at T_ref), `reference_temperature` (default the
  `[fluid] temperature`) and `reference_density`. The solved pressure is
  the reduced pressure p′ = p − ρ_ref g·(x − x_ref); the output has the
  static pressure p as `static_pressure` as well. The values of the
  pressure boundaries are p′ by default; `pressure_kind = "static"` reads
  them as static pressures and converts them per face. Works with SIMPLEC
  and the coupled algorithm (the kernels are shared by CPU and GPU), with
  the energy equation and with solid regions (the fluid temperature drives
  the force).
- **Balanced force (D154).** The buoyancy enters every face as the
  pressure gradient does: the same face differences, the same
  least-squares operator in the cells, the matching Rhie–Chow term and
  boundary pressure. A closed box with a stable stratification along an
  oblique gravity stays at rest to 10⁻¹⁵ m/s on perturbed hex, tet, poly
  and mixed meshes, with SIMPLEC and the coupled algorithm (V7.7 (a)).
- **Buoyancy pseudo-time step.** Steady runs with buoyancy take a
  pseudo-time step 1/N (N the buoyancy frequency, scaled by
  `[solver.pseudo_time] factor`) in SIMPLEC and SIMPLE, and the coupled
  algorithm's automatic step takes the buoyancy velocity as a lower bound.
  Without it a stable stratification on a coarse mesh never settled.
- **Kinetic-energy ledger** counts the work of the buoyancy force.
- **Codes** BUO001 (no gravity), BUO002 (Boussinesq out of range), BUO003
  (no energy equation or no β), BUO004 (gravity with energy but no model),
  BUO006 (a static-kind pressure value over a tall opening); the
  variable-density model is refused with ENE006 until stage 7c.
- **Cases and docs.** `buoyancy_enclosure_plume`, `buoyancy_box_at_rest`,
  `buoyancy_slot` (analytic `buoyant-slot`, observed orders 2.14 and
  2.13), `buoyancy_channel_reduced`, `buoyancy_channel_static`,
  `buoyancy_open_box`; gallery and recipe 186 ("Switch on gravity and
  buoyancy; which pressure do my boundary values mean").
