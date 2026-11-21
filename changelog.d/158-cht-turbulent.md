- **Turbulent conjugate heat transfer (stage 7d5).** Solid regions now
  work with every turbulence model. The fluid side of each fluid–solid
  interface is the thermal wall function of the flow (Kader or
  Jayatilleke, D135), h_w = ρ c_p u*/T⁺, in series with the solid's half
  cell, q = (T_P − T_N)/(1/h_w + δ_N/k_s); y* and μ of the first fluid
  cell are taken in every outer iteration, the fluid conducts with
  k + c_p μ_t/Pr_t, and μ(T) goes back to the laminar viscosity of the
  turbulence model. Check V7.36: a turbulent channel over a conducting
  wall gives the same heat flow as the fluid-only path with the wall
  function in series with a convection condition, to 10⁻⁴ for k-ε and
  SST (`cht_turbulent_channel.toml`). Gallery and recipe 158: a finned
  aluminium heat sink in a turbulent air channel
  (`cht_heat_sink_turbulent.toml`).
- **Fixes.** A `[[schedule]]` entry no longer drops the energy and scalar
  sub-tables of the boundaries (an energy case with an inlet and a
  schedule was refused). A temperature that turns NaN in a diverging run
  is reported as a divergence instead of stopping in the property
  interpolation.
