# Engineering gallery cases (`examples/applications/`)

Part of the [Case catalogue](../cases.md).

| Case | Mesh | Shows | Where |
|---|---|---|---|
| `applications/stirred_tank/case.toml` | Gmsh tetrahedra from `mesh.py --level 1`, 607 032 cells | the 1 m tank with a three-blade hydrofoil, three baffles, sparger and feed pipe, moving reference frame, k-ω SST, 90 rpm; `post.py` gives the power and pumping numbers, `images.py` the figures | [application 131](../applications/rotating-machinery.md#stirred-tank) |
| `applications/stirred_tank_fyfe/case.toml` | Gmsh tetrahedra from `mesh.py --level 0`, 194 827 cells | the 2 L lab tank of Fyfe et al. (four-blade pitched-blade turbine, four baffles) at 200 rpm, steady frozen rotor: the power number against measured values; `geometry_figure.py` draws the section | [application 161](../applications/rotating-machinery.md#fyfe-stirred-tank) |
| `applications/stirred_tank_fyfe/tracer.toml` | the same mesh | the transient tracer run of the tank, restarted from the steady flow: a 20 mL blob at the liquid surface, k-ω SST, Sc 1 000, Sc_t 0.7; `mixing_time.py` gives θ₉₅ against the measured plateau time of the paper and the Grenville correlation | [application 161](../applications/rotating-machinery.md#fyfe-stirred-tank) |
| `applications/shell_tube_hx/mesh.toml` | Gmsh tetrahedra from `mesh.py --level 2`, 983 249 cells | the shell and tube heat exchanger: six cell zones (hot air, water, tubes, tubesheets, baffles, housing), four nozzle boundaries, conformal interfaces; `foehn mesh` checks it, `images.py` draws the section and the cut-away | [application 203](../applications/heat-exchangers.md#shell-tube-hx) |
| `applications/shell_tube_hx/case_coflow.toml`, `case_counterflow.toml` | the same mesh (`mesh.py --level 2`) | the heat exchanger in steady conjugate heat transfer: hot air (300 °C, 15 m/s, 1.5 bar) on the shell side, water (20 °C, 5 m/s, 1 bar) in coflow or counter flow on the tube side, carbon steel, housing at 25 °C, k-ω SST; waits for several fluid materials in one case (D194) | [application 203](../applications/heat-exchangers.md#shell-tube-hx) |
