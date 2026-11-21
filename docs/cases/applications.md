# Engineering gallery cases (`examples/applications/`)

Part of the [Case catalogue](../cases.md).

| Case | Mesh | Shows | Where |
|---|---|---|---|
| `applications/stirred_tank/case.toml` | Gmsh tetrahedra from `mesh.py --level 1`, 607 032 cells | the 1 m tank with a three-blade hydrofoil, three baffles, sparger and feed pipe, moving reference frame, k-ω SST, 90 rpm; `post.py` gives the power and pumping numbers, `images.py` the figures | [application 131](../applications/rotating-machinery.md#stirred-tank) |
| `applications/stirred_tank_fyfe/case.toml` | Gmsh tetrahedra from `mesh.py --level 0`, 196 713 cells | the 2 L lab tank of Fyfe et al. (four-blade pitched-blade turbine, four baffles) at 200 rpm; `compare.py` draws the planes and profiles against the paper's published CFD plane | [application 161](../applications/rotating-machinery.md#fyfe-stirred-tank) |
