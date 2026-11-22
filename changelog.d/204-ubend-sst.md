- **U-bend (D45).** New cases `ubend_sst.toml` and `ubend_sst_sharp.toml`:
  the two-dimensional U-duct of Monson et al. (1990), Rc/H = 1, at
  Re_b = 10⁵ with k-ω SST and a developed channel profile at the inlet (new
  case `ubend_channel_profile.toml`, table by
  `examples/studies/ubend_inlet_profile.py`), on structured meshes from
  `examples/meshes/ubend.py` with the rounded turn of the experiment or with
  sharp edges (`--sharp`: a square divider end and end wall). The rounded
  case converges on three meshes from a potential-flow start (GCI in
  gallery 204, `examples/studies/ubend_gci.py`): the inner wall separates at
  θ = 139.8° (GCI 0.19 %) and reattaches 2.19 H behind the turn (GCI 2.9 %);
  the wall pressures have a GCI of 0.27 % or less. The curvature correction
  of SST moves the separation to 127.5° and the reattachment to 2.57 H. The
  sharp-edged duct converges on the two coarser meshes with a gentler start
  (relaxation 0.5, linear upwind from iteration 601, relaxation 0.8 from
  1 201). No comparison with measurements yet (pending). Gallery and
  recipe 204.
