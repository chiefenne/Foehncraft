- **Asymmetric plane diffuser (D45).** New cases `diffuser_sst.toml` and
  `diffuser_sst_sharp.toml`: the diffuser of Buice & Eaton (10°, area ratio
  4.7) at Re_H 20 000 (centreline velocity) with k-ω SST and a developed
  channel profile at the inlet (new case `diffuser_channel_profile.toml`,
  table by `examples/studies/diffuser_inlet_profile.py`), on structured
  meshes from `examples/meshes/diffuser.py` with the corners rounded as in
  the experiment (4.3 H) or sharp (`--sharp`, cells shortened towards the
  corners). Both converge on three meshes (GCI in gallery 199,
  `examples/studies/diffuser_gci.py`): with rounded corners the reattachment point (28.94 H on 111 360 cells) has a
  GCI of 0.36 %, the separation point (2.20 H) 4 % and the pressure recovery
  0.03 % or less. No comparison with
  measurements yet (pending). Gallery and recipe 199.
- **Why the diffuser did not converge before.** A sharp wall corner makes
  the least-squares gradient of thin wall cells wrong by O(1), and linear
  upwind from the first iteration diverged on it; trapezoidal wall cells
  switched the D50 remediation on, on the coarse levels only; a uniform
  start broke continuity in the widening; and relaxation 0.5 only looked
  like a stall. The case now starts from potential flow with upwind for
  200 iterations, runs at relaxation 0.8, and judges convergence on probes
  in the log (`docs/pitfalls.md`).
- **The momentum scheme may change during a run.** `solver.momentum.scheme`
  may be scheduled (`[[schedule]]`) or re-read (`control.toml`) between
  `"upwind"` and `"linear-upwind"`: a first-order start, then second
  order. A change to or from `"central"` is refused.
- **A transient run may start from a steady restart file.** `--restart`
  with a steady restart file and a transient case starts the transient
  run from the steady fields at `time.start_time` (no older time levels);
  the log says so.
