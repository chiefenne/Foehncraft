# Heat exchangers

Part of the [Engineering Gallery](../applications.md): shell and tube,
plate and finned exchangers, with the fluids and the metal in one mesh for
conjugate heat transfer (recipes [156](../recipes/energy.md#cht-solids) to
[190](../recipes/energy.md#cht-cycle)).

| # | Application | Sector | Confidence |
|---|---|---|---|
| 203 | [Shell and tube heat exchanger: a conformal multi-region tet mesh from a STEP assembly](#shell-tube-hx) | heat exchangers | mesh only |

---

<a id="shell-tube-hx"></a>
## 203. Shell and tube heat exchanger: a conformal multi-region tet mesh from a STEP assembly

**Problem.** A compact shell and tube exchanger: hot air flows around the
tube bundle on the shell side, guided by three segmental baffles, and heats
water that flows through seven tubes from one end chamber to the other. To
compute the heat duty both fluids and every metal part must share one mesh
whose faces match across every interface. This entry builds that mesh from
the CAD assembly with Gmsh; the flow run follows with the boundary
conditions.

| | |
|---|---|
| Sector | heat exchangers |
| Geometry | `examples/applications/shell_tube_hx/shell_tube_hx.step` (realistic, not from a real project): shell of inner diameter 90 mm, wall 6 mm, dished heads, 474 mm long; seven tubes 21 × 3 mm, 300 mm long; two tubesheets and three segmental baffles 6 mm thick (cut 20 mm off the axis, alternating); tube-side nozzles DN 40 on top, shell-side nozzles DN 30 on the sides |
| Mesh | Gmsh tetrahedra (HXT, Netgen optimised), no prism layers; level 2 (3 mm): 983 249 cells, 196 435 nodes; largest non-orthogonality 73.6° (2 faces above 70°), 74 faces of skewness above 0.9 (D50 remediation); level 0 (5.1 mm): 477 921 cells |
| Cell zones | `shell_side` (hot air), `tube_side` (water: the end chambers and the tube bores), `tubes`, `tubesheets`, `baffles`, `housing` |
| Boundaries | `tube_nozzle_front`, `tube_nozzle_rear` (water), `shell_nozzle_front`, `shell_nozzle_rear` (air), `housing_outer` (every outer face of the housing) |
| Interfaces | found by the solver from the cell zones (`interface-<region>`, recipe 156); the two fluids touch nowhere |
| Run cost | meshing 40 s (level 0) to 2 min (level 2) on one core |
| Confidence | **mesh only**: C6 passes (`foehn mesh`: valid, the boundary encloses the volume to 2·10⁻¹⁶, no cell left unremediated); the flow checks follow with the run |
| Case files | `applications/shell_tube_hx/mesh.toml` (mesh check); the flow case follows |

### The mesh

<img src="../figures/applications/hx_mesh_section.png" alt="Section through the mesh of the shell and tube heat exchanger: hot air on the shell side, water in the end chambers and the tube bores, tubes, tubesheets, baffles and housing in their own colours" width="800" height="400">

*The section x < 0 (whole cells): water in the end chambers and in the
bore of the centre tube (blue), hot air between the tubesheets (orange),
the three baffles alternating top and bottom, the housing around them.*

<img src="../figures/applications/hx_mesh_cutaway.png" alt="Cut-away of the housing with the surface mesh of the tube bundle, tubesheets and segmental baffles" width="800" height="400">

*The housing cut open, the bundle whole: seven tubes, the two tubesheets
and the three segmental baffles with their surface mesh.*

### Volumes and interfaces

The zone volumes match the CAD solids: tubes 7 × 50 893.8 mm³, tubesheets
2 × 23 623.2 mm³, baffles 3 × 18 470.4 mm³ exactly, the housing to
1·10⁻⁵. The interfaces the solver will couple:

| Zone | Zone | Area (cm²) |
|---|---|---:|
| shell_side | tubes | 1 263.4 |
| shell_side | housing | 883.8 |
| shell_side | baffles | 193.0 |
| shell_side | tubesheets | 78.7 |
| tube_side | tubes | 1 013.4 |
| tube_side | housing | 643.5 |
| tube_side | tubesheets | 78.7 |
| tubes | baffles | 66.7 |
| tubes | tubesheets | 55.4 |
| tubesheets | housing | 33.9 |
| baffles | housing | 32.9 |

### How the mesh is built

`mesh.py` reads only the housing from the STEP file. The tubes, tubesheets
and baffles are rebuilt as exact OpenCASCADE primitives with the
dimensions of the STEP file (the same to 10⁻¹² mm): a Boolean cut of the
STEP tubes from the void inside the housing returned wrong pieces (tube
walls classified as fluid), the primitives do not. One `fragment` of a box,
the housing, four nozzle plugs and every internal part then gives each
piece once, with shared faces. A piece takes its zone from the inputs that
contain it (housing before tube bore before tube wall before tubesheet
before baffle); a piece only inside the box is fluid, on the shell side
between the tubesheets and on the tube side beyond them. The plugs, short
discs just outside the nozzle ends, close the nozzles: their flat faces
become the four inlet and outlet boundaries. The script stops when the two
fluids share a face or a nozzle face is missing.

### Limits

The flow run needs two fluids with their own materials in one conjugate
case (air on the shell side, water on the tube side); this version takes
one `[fluid]` material for the energy equation and refuses a second
(ENE006), so the run waits for that solver step. No prism layers yet: the
first cell at the tube walls is about 1 mm, so a turbulent shell-side run
uses wall functions. One mesh family (levels 0 to 3, ratio 1.3) is ready
for the GCI.

Recipe: [203](../recipes/applications.md#shell-tube-hx-recipe).
