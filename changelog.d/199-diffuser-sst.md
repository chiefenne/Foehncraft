- **Asymmetric plane diffuser (D45).** New case `diffuser_sst.toml`: the
  diffuser of Buice & Eaton (10°, area ratio 4.7) at Re_H 20 000 (centreline velocity) with k-ω SST and a developed
  channel profile at the inlet (new case `diffuser_channel_profile.toml`,
  table by `examples/studies/diffuser_inlet_profile.py`)
  on a structured mesh from `examples/meshes/diffuser.py`. The coarse mesh
  (5 280 cells) gives a steady solution with separation at 2.3 H and
  reattachment at 29.3 H; the finer meshes (21 120 and 84 480 cells) do not
  reach a steady state with the steady solvers, so there is no grid-convergence
  result, written as a gap. No comparison with measurements yet (pending).
  Gallery and recipe 199.
