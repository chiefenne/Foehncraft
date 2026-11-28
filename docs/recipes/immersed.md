# Immersed bodies

Part of the [Recipes](../recipes.md): put bodies into the flow without meshing them (L13). The keys are on the [key page](../keys.md) (`immersed.*`).

| # | Recipe |
|---|---|
| 270 | [Put a body in the flow without meshing it](#penalised-body) |
| 291 | [Search the shape of a duct with a density field](#topology-search) |
| 313 | [Cut a body out of a mesh](#cut-body) |
| 310 | [Check an STL surface before it becomes a body](#surface-check) |
| 315 | [Cut a mesh exactly by an STL surface](#cut-moments) |
| 322 | [Put an STL body into the flow](#stl-body) |
| 357 | [Run a turbulent flow with wall functions on a cut body](#cut-wall-functions) |
| 359 | [Keep a cut body as a solid region (conjugate heat transfer)](#cut-cht) |
| 360 | [Give a cut body inlets, outlets and a turning wall](#cut-patches) |
| FFF | [Refine the mesh near a body](#refine) |

---

<a id="penalised-body"></a>
## 270. Put a body in the flow without meshing it

**Add an `[immersed]` table with one `[[immersed.body]]` per body** to a laminar, isothermal flow case on any mesh. The mesh needs no boundary for the body; the body may cross cells and faces anywhere:

```toml
[immersed]
permeability = 1e-8        # K (m²): sink μ χ V / K; no default
# samples = 4              # sample points per direction and tetrahedron (1 to 16)

[[immersed.body]]
name = "cylinder"
shape = "cylinder"         # "cylinder", "sphere" or "box"
centre = [0.2, 0.2, 0.0]   # a point on the axis
radius = 0.05
# axis = [0.0, 0.0, 1.0]   # default z
# length = 0.1             # default endless

[[immersed.body]]
name = "block"
shape = "box"
lower = [0.6, 0.1, -1.0]
upper = [0.7, 0.2, 1.0]

[output]
derived = ["solid-fraction"]   # χ in fields.vtkhdf
```

**Choose K.** Smaller K is closer to a no-slip wall and stiffer, but it does not make the wall sharper: the solid fraction of the cells along the wall smears the body outward by up to half a cell, and a large K lets the flow seep in over the depth √K. Take K about 0.01 h² (h the cell width at the body): on the sweeps of [gallery 291](../gallery/immersed.md#topology-bend) the drag and mass flow errors are then a few per cent and fall with h at about first order; 10⁻³ h² and below is two to three times as far off on the sphere (the cylinder case uses 4·10⁻⁴ h²). SIMPLEC takes any K; the coupled solver converges from rest with a sink up to about 10⁵ times the viscous diagonal (μ/K below about 10⁵ μ/h²), above that its linear solver can break down. Check the effect of K by halving it once.

**Run and read the forces.**

```sh
target/release/foehn flow examples/cases/flow/immersed_cylinder_re20.toml --out output/imm
```

The log lists every body with its solid volume on the mesh against the exact volume (IMM104; IMM103 warns when the mesh does not resolve the body). `result.json` has, per solve, `immersed`: name, solid volume, exact volume, cells and the force of the fluid on the body (N). `[verification] benchmark = "schafer-turek-2d1"` takes the body named `cylinder` (or `benchmark_boundary`) when the mesh has no boundary of that name.

**Limits.** Laminar and low-Re only, and isothermal: a turbulence model, the energy equation, passive scalars, gravity, regions, rotating zones and sprays are refused (IMM102). The wall is smeared over about a cell: drag and pressure come out well, wall shear stress and lift do not ([gallery 270](../gallery/immersed.md#penalised-cylinder)). Use it for a first look and topology; a sharp wall comes with the cut cells of L13 stage 2.

---

<a id="topology-search"></a>
## 291. Search the shape of a duct with a density field

**Add `[immersed.topology]`** to a laminar, steady, isothermal case: the cells whose centroid lies in the box `lower`–`upper` become design cells with a fluid fraction γ (1 fluid, 0 solid); the search minimises the power the flow takes from its boundaries, with `fluid_fraction` of the box's volume fluid. Bodies of `[[immersed.body]]` stay fixed solid beside the design.

```toml
[flow]
convection = false           # Stokes flow: the gradient is exact

[immersed]
permeability = 4e-5          # K of a fully solid cell: μ/K = 2.5e4 here

[immersed.topology]
lower = [0.0, 0.0, -1.0]     # the design region
upper = [1.0, 1.0, 1.0]
fluid_fraction = 0.25
# interpolation = 0.1        # q of χ = q (1 − γ)/(q + γ); smaller: grey cells penalised harder
# filter_radius = 0.0        # m; 1.5 h against checkerboards and mesh dependence
# move_limit = 0.2           # largest change of γ per update
# every = 20                 # flow iterations between two updates
# designs = 100              # most updates
# change = 0.01              # stop when no γ changes by more
# start = 0.25               # start value of γ (default fluid_fraction)

[output]
derived = ["solid-fraction"] # 1 − γ̃ in fields.vtkhdf
```

**Choose μ/K** large against the viscous term (μ/K ≥ 10⁴ μ/L² for the domain size L), so that solid cells carry no flow; the coupled solver needs the sink below about 10⁵ times its viscous diagonal from rest, SIMPLEC takes any K.

**Run it** like any case:

```sh
target/release/foehn flow examples/cases/flow/immersed_topology_bend.toml --out output/bend
```

The log has one `topology design` line per update (Φ, fluid fraction, largest change, grey level); `topology.csv` the same rows; `topology.stl` the surface γ = ½ (marching tetrahedra, closed and manifold), ready for a body-fitted mesh or the cut cells of stage 2; `result.json` per solve `topology` (designs, converged, first and final Φ, fluid fraction, grey level, triangles, surface area) and the force on the design as the body `design`.

**What the gradient is.** ∂Φ/∂γ = (μ/K) χ′(γ) |u|² V, the derivative of the dissipation of Stokes flow (self-adjoint): exact for `convection = false`, an approximation with convection (0.5 % at Re 20, larger at higher Re; the adjoint solver L17 is not built yet). Checked against central differences in `tests/immersed.rs`.

**Limits.** Steady, laminar and isothermal only (IMM102 refuses a transient run and the settings of recipe 270), one rank, a box as design region. The result is a topology, not a final shape: the wall is smeared over a cell, so check the design on a body-fitted mesh made from `topology.stl`.

---

<a id="cut-body"></a>
## 313. Cut a body out of a mesh

**Set `method = "cut"`** in `[immersed]` and name a wall condition for every body: the cut takes the bodies out of the mesh before the run, and each body surface becomes a wall boundary with the body's name.

```toml
[immersed]
method = "cut"

[immersed.cut]
min_volume_fraction = 0.3   # merge cut cells below this fraction (0 to 0.9)
# edge_root = "exact"       # crossing points on the true surface; "linear": marching cubes
# wall_point = true         # wall faces fanned around a point on the surface
# face_point = true         # cuts across the faces bent onto the surface

[[immersed.body]]
name = "cylinder"
shape = "cylinder"          # "cylinder", "sphere", "box", "field" or "stl"
centre = [0.2, 0.2, 0.0]
radius = 0.05

[boundary.cylinder]
type = "wall"               # any wall condition: velocity, temperature, heat flux
```

No permeability is needed. The log has the cut report (IMM203): cells cut, removed, split and merged, the smallest kept volume fraction, the closure of the cut cells and the wall area of each body against the exact area where it is simple. IMM204 warns when a small cell found no neighbour to merge into or a face is left that the discretisation refuses; refine the mesh near the body or move it by a small fraction of a cell. Forces and coefficients come from the body boundary like from any wall (`boundaries.csv`, `[verification] benchmark_boundary`).

**Choose the merge threshold.** On the Re 20 cylinder the lift moves from 4 % to 0.4 % below the benchmark range between 0.1 and 0.4, drag and Δp hardly move ([gallery 313](../gallery/immersed.md#cut-cylinder)); 0.3 is the default. Merged cells are larger and non-convex, so do not go far above 0.4.

**Resolution.** The body is seen through its values at the mesh nodes: a part thinner than a cell can vanish, and the wall faces are flat between crossing points (second order in the cell width). Put about 20 cells across the smallest feature.

**Limits.** STL bodies: [recipe 322](#stl-body). Regions and conjugate heat transfer, a topology search and rotating zones are refused (IMM202). One rank only.
<a id="surface-check"></a>
## 310. Check an STL surface before it becomes a body

**Run the check** on the file (ASCII or binary STL; `scale` converts the file's unit to metres):

```sh
echo '{"input": "body.stl", "scale": 1e-3}' | target/release/foehn call mesh.surface --input -
```

It prints the counts (triangles, distinct corners, pieces, area, enclosed volume, closed or open) and one line per finding; with `--json` the whole report (`clean`, the counts and the first 20 defective triangles or pairs of each kind). Nothing is repaired.

| Code | What | What to do |
|---|---|---|
| IMM001 | no triangles, non-finite corners | export again |
| IMM002 | holes (edges with one triangle) | small holes are fine for inside/outside (the winding number closes them); close large ones in CAD |
| IMM003 | edges with more than two triangles | separate touching bodies |
| IMM004 | inconsistent orientation | fix the normals in CAD; a few flipped triangles only dent the inside test locally |
| IMM005 | self-intersections | tolerated by the winding number; a cut (stage 2) needs a clean surface |
| IMM006, IMM007 | degenerate, duplicate triangles | remove in CAD |
| IMM008 | inward normals (negative volume) | flip all triangles |
| IMM009 | box larger than 20 m, near-duplicate corners | set `scale`; weld in CAD |

**In code** (Rust, `foehncraft_mesh::surface`): `Surface::read_stl(path, scale)`, `check(&s)`, `s.winding_number(x)` (fast, β = 2) or `winding_number_exact(x)`, `s.inside(x)`, `s.signed_distance(x)` and `signed_distances(&points, Some(band))` (negative inside, ±band beyond the band), `s.crossings(p, q)` and `first_crossing(p, q)` (the exact points where a segment enters or leaves the body: the root along a mesh edge for the cut of stage 2).

**Limits.** The fast winding number is within about 4e-3 of the exact sum next to the surface, which decides inside and outside except within a hair of the surface (there `inside` re-evaluates more precisely). The signed distance of a self-overlapping surface counts the buried triangles (its sign is right, its magnitude inside is the distance to the nearest triangle). The signed-distance benchmark and the exact cut followed in [recipe 315](#cut-moments).

<a id="cut-moments"></a>
## 315. Cut a mesh exactly by an STL surface

The geometry kernel gives the numbers a cut-cell method needs, for any polyhedral mesh (Rust, `foehncraft_mesh::surface`):

```rust
let s = Surface::read_stl(path, 1e-3)?;            // file in mm
let geo = Geometry::compute(&mesh);
let cut = s.cut_moments(&mesh, &geo);
// cut.cell_inside_volume[c], cut.cell_inside_centroid[c]   part of cell c inside the body
// cut.face_inside_area[f], cut.face_inside_centroid[f]     part of face f inside (aperture = area / |S_f|)
// cut.walls: WallPiece { cell, triangle, area (out of the body), centroid }
```

The fluid part of a cell is its volume minus the inside volume; the fluid aperture of a face is one minus the inside one. Use the signed distance (recipe 310, `signed_distances(&nodes, Some(band))`) where a level set is enough, and `cut_moments` where volumes, apertures and wall faces must be exact: the cut of a sphere is 40 to 2500 times more accurate than a linear level-set cut on 64³ to 8³ cells (gallery 315), and conservative (both cells of a face see the same aperture, every cell closes).

**Which signed distance.** The default (closest point with winding-number sign) won the benchmark of gallery 315: exact on clean surfaces and sharp edges, no sign errors on holes and soups. Check the surface first (recipe 310): a patch of flipped triangles (IMM004) gives wrong signs next to it.

**Limits.** The surface must be in general position with respect to the mesh except for edges lying in a mesh face (handled); a surface lying flat in a mesh face or a surface corner exactly on a mesh node is not decided consistently yet. The cut gives moments, not polygons: building cut cells into a mesh is stage 2 (D227).
**Limits.** The fast winding number is within about 4e-3 of the exact sum next to the surface, which decides inside and outside except within a hair of the surface (there `inside` re-evaluates more precisely). The signed distance of a self-overlapping surface counts the buried triangles (its sign is right, its magnitude inside is the distance to the nearest triangle). Not yet: the signed-distance benchmark against fast sweeping and the heat method, and the exact triangle–cell clipping (the next parts of stage 0).

<a id="stl-body"></a>
## 322. Put an STL body into the flow

**Name the file** in a body of `shape = "stl"`; the surface becomes a wall boundary with the body's name when it is cut out of the mesh, or a penalised body without `method = "cut"`:

```toml
[immersed]
method = "cut"              # or leave it out (penalisation, needs permeability)

[[immersed.body]]
name = "part"
shape = "stl"
file = "part.stl"           # ASCII or binary, relative to the case file
scale = 1e-3                # file in millimetres (default 1)

[boundary.part]
type = "wall"
```

**Check the surface first** ([recipe 310](#surface-check)): the run does the same check and writes its findings to the log (IMM002 to IMM009); holes and overlaps are tolerated through the winding number, a patch of flipped triangles gives wrong signs next to it. A file that cannot be read (missing, no triangles, a corner that is not a number) stops the run with CASE006.

**Sharp edges and corners stay sharp by default.** Where the triangles under a wall loop meet at more than 15°, the wall faces are fanned around the point of dual contouring on the edge or corner, and `face_point = true` bends the cuts across the faces onto the surface; on a box turned against the mesh this makes the wall area seven to ten times closer at 32³ than the closest point ([gallery 322](../gallery/immersed.md#stl-cut)). Where such a point would give a face the discretisation refuses, the cut drops it for the cells around it (IMM203 says how many cells, in how many rounds); `face_point = false` and `wall_point = false` switch them off everywhere.

**Resolution.** About 6 cells per radius give the drag of a sphere array within 2.3 %, 10 cells within 0.9 % (hexahedra, second order). A part of the surface thinner than a cell can vanish; put two or more cells across thin walls.

**Limits.** Hexahedral meshes are the best checked; on tetrahedra the force balance of the wall is within 5 %, on polyhedra the cut is less accurate (the faces are not planar). Merged cells keep faces that do not lie in one plane apart, so two cells can share several faces (IMM203 counts them). Wall functions: [recipe 357](#cut-wall-functions), conjugate heat transfer: [recipe 359](#cut-cht). Other wall conditions per surface patch and near-body refinement follow (stages 2e and 2f). One rank only.

---

<a id="cut-wall-functions"></a>
## 357. Run a turbulent flow with wall functions on a cut body

**Nothing extra to set.** A cut body (`method = "cut"`) in a turbulent case with wall functions (any k-ε treatment but the low-Re ones, SST with its automatic wall treatment) gets them on its wall: every cut or merged wall cell takes the exact distance of its centroid to the body, and its wall function takes the velocity from the cell one layer out along the wall normal.

```toml
[immersed]
method = "cut"

[immersed.cut]
# wall_sample = true        # sample one layer out (default); false: the cut cell itself

[[immersed.body]]
name = "plate"
shape = "plane"             # a half-space: an inclined wall across the mesh
centre = [0.0, 0.02, 0.0]   # a point on the plane
normal = [-0.0349, 0.9994, 0.0]   # out of the solid into the fluid

[boundary.plate]
type = "wall"

[turbulence]
model = "k-epsilon"
wall_treatment = "scalable"
```

**Read the y⁺ of the body.** LIM011 in the log and in `result.json` gives the smallest, mean and largest first-cell y⁺ of every cut body; `[output] wall_csv = true` writes the `y_plus` of every wall face (`wall_<body>.csv`), the y⁺ map. The y⁺ of a cut wall spreads more than on a fitted wall (a plate crossing cell rows: 30 to 140 where the fitted plate has 100), because the cut cells have every height between a fifth and one and a half cells.

**Mesh for the sample.** The sampled cell has to lie in the log layer: keep the first two or three cell rows near the body below y/δ ≈ 0.15. With 2 mm cells on a plate with a 16 mm boundary layer the skin friction came out 17 % high; with 1 mm cells within 7 % ([gallery 357](../gallery/immersed.md#cut-plate)). Where the mesh is too coarse for that, `wall_sample = false` takes the cut cell itself: the mean is then right but single faces scatter by ±40 % or more.

**Start the flow along the wall.** A uniform start velocity that points into a cut body at high Reynolds numbers can diverge in the first iteration (net inflow into the cut cells): start with the velocity along the body or from a coarser run.

**Limits.** Resolved walls (low-Re treatments, ζ-f to the wall) take the exact distance but no sample. Single wall faces scatter about ±15 % in C_f even with sampling; integrated forces and the mean along the body are closer. One rank only.

---

<a id="cut-cht"></a>
## 359. Keep a cut body as a solid region (conjugate heat transfer)

**Name the body as the zone of a solid region.** The cut then keeps the inside of the body as the cells of that region; the cut faces become the region interface (D160), the flow sees the wall `interface-<region>`:

```toml
[immersed]
method = "cut"

[[immersed.body]]
name = "pin"
shape = "cylinder"          # any shape, an STL body too
centre = [0.2, 0.2, 0.0]
radius = 0.05

[materials.aluminium]
kind = "solid"
density = 2700.0
specific_heat = 900.0
conductivity = 200.0

[regions.pin]
kind = "solid"
zones = ["pin"]             # the body, not a mesh zone
material = "aluminium"

[energy]
initial = 300.0

[[energy.sources]]
region = "pin"
power_density = 2e5         # W/m³
```

The body needs no `[boundary.<body>]`: its wall is the interface. A body without a region is cut out as before (a wall); both kinds can be in one case. The energy equation is required (IMM202 otherwise).

**Read the cut.** The log's mesh line gives, per solid body, the interface area and the solid volume against the exact ones (`cut body pin (solid region pin): interface area …, solid volume …`) and the solid cells, the small ones merged and the smallest kept fraction (`cut cells both ways`). The heat balance of the energy line closes to round-off when the energy solve has converged.

**Converge the energy equation.** The solid with k_s ≫ k_f is stiff in one matrix with the fluid; give the energy solve its own tolerance and enough linear iterations:

```toml
[solver.energy]
relaxation = 1.0
tolerance = 1e-10
max_iter = 300
```

Without it the pin of [gallery 359](../gallery/immersed.md#cut-cht) stalled with a heat imbalance of 33 %.

**Two cells across, or a thin wall.** A solid body must be about two cells thick where it is thinnest; a thinner one is refused with IMM205. Refine the mesh in the body, or model the wall as a thin wall on a region interface (`[interfaces.<name>] thin_wall`, D189), which carries its conduction without cells.

**Limits.** Hexahedra and tetrahedra are the best checked; on polyhedra a few flat solid slivers stay. Solids that move (stage 5) and the thin-wall model on a cut interface are not built. One rank only.

---

<a id="cut-patches"></a>
## 360. Give a cut body inlets, outlets and a turning wall

**Name the patches in the STL file.** Write the surface as an ASCII STL file with one `solid <name>` per patch (most CAD tools export the faces of a part as named solids), or tag the triangles of a binary file with their 16-bit attribute. Then:

```toml
[immersed]
method = "cut"

[[immersed.body]]
name = "cylinder"
shape = "stl"
file = "blown_cylinder.stl"
patches = "solids"              # ASCII: one patch per solid
# patches = { slot = 1, wall = 2 }   # binary: patch = attribute

[boundary.cylinder-slot]        # <body>-<patch>
type = "velocity-inlet"
normal_velocity = 0.4
[boundary.cylinder-wall]
type = "wall"
```

Any condition goes on a patch: inlets, outlets, walls with their own wall velocity, roughness or energy condition. Faces of triangles without a listed tag stay in the boundary `<body>`. The log lists the faces and area of every patch.

**Turn the wall.** `omega` (rad/s) turns the wall of a body about its `axis` through its `centre` (a box: its middle):

```toml
[[immersed.body]]
name = "rotor"
shape = "cylinder"
centre = [0.0, 0.0, 0.0]
radius = 0.05
omega = 10.0                    # the wall moves with ω × (x − centre)

[boundary.rotor]
type = "wall"                   # no `velocity`: omega gives it
```

The boundary summary of every wall of the body (`result.json`, the log) then has its `torque` about the axis; give `omega = 0` for the torque of a body at rest. The body itself does not move: this is exact for a body of revolution about its axis (a shaft, a disc, a cylinder) and for a body at rest. A wall `velocity` of your own and `omega` on the same body are refused (IMM101).

**Keep the fluid inside.** `invert = true` turns a cylinder, sphere or box inside out: the fluid inside, the solid outside. A pipe or the outer wall of a gap then needs no fitted mesh. Boundaries of the mesh that lie wholly in the solid are left out and need no condition (a condition the case gives them anyway is kept, on an empty boundary).

**Limits.** A wall face of a cut cell is not split along a patch edge: patch edges are resolved to the cell size. Solid regions (recipe 359) cannot turn or have patches. One rank only.
**Limits.** A wall face of a cut cell is not split along a patch edge: patch edges are resolved to the cell size. Solid regions (recipe 359) cannot turn or have patches. One rank only.

---

<a id="refine"></a>
## FFF. Refine the mesh near a body

**Start from a coarse box mesh and let the bodies refine it.** Add `[immersed.refine]` to a case with `[[immersed.body]]` bodies (cut or penalised):

```toml
[mesh.generate]
type = "hex"
dimension = 2
upper = [2.2, 0.41, 0.01]
cells = [220, 41, 1]           # h = D/10 away from the body

[immersed]
method = "cut"

[immersed.refine]
levels = 2                     # halve the cells near the bodies twice: D/40
distance = 0.08                # the finest level within 0.08 m of a surface,
                               # the next one within 0.16 m, ...
# layers = 4                   # or: each level 4 of its own cells wide (default)

[[immersed.body]]
name = "cylinder"
shape = "cylinder"
centre = [0.2, 0.2, 0.0]
radius = 0.05
```

The hexahedra of the mesh near the surfaces are halved before the cut, 2:1 balanced, into hexahedra with hanging-node faces; a planar mesh (one cell thick) is refined in its plane. The log (IMM206) gives the cells per level, the faces between levels and the largest non-orthogonality and skewness.

**Choose the widths.** Make the finest band wide enough for the flow features at the body: the lift of the Schäfer–Turek cylinder ([gallery FFF](../gallery/immersed.md#refine)) is 5 % below the reference range with the finest band 0.02 m (0.2 D) wide and inside it at 0.08 m (0.8 D); drag and pressure difference were within 0.2 % of their ranges already at 0.02 m. Use `distance` for a convergence study (`[verification] levels`): the refined region is then the same in metres on every level, while `layers` narrows it with the cells.

**Limits.** Only hexahedra are refined; tetrahedra, prisms and polyhedra are kept whole. The refinement follows the distance to the bodies only, not the solution, and never coarsens. A mesh with periodic pairs is refined before the pairing. One rank only.
