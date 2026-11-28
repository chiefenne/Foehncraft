# L13 – Immersed boundaries: staged plan

Plan only, no code (maintainer choice "Plan first", 2026-11-22). Builds on
[step 26](step26-design-optimization.md) (§3 methods, §8 layered idea, §9 2D
prototype), the amended L13 and D58, D73, D108, D211 in
[requirements.md](../requirements.md). Clean room (D201): published
literature only, open-source code may be read where the licence allows, never
copied, no commercial solver names. Sizes are agent time: code in minutes,
runs in hours, CI separate; "not sized yet" where unsure.

Ids: D225 to D234 (next free after the adjoint plan was D225 and L52),
entered in [requirements.md](../requirements.md) with this page; the table in
§0 maps them to the stages.

**Decisions of 2026-11-22.** The maintainer answered the eleven questions of
§9 (answers in the table there, one line each): no accuracy numbers yet (Q1),
tolerant dirty STL with Crane's methods as candidates (Q2), merging of small
cells (Q3), thin solids refused (Q4), penalised walls laminar and low-Re only
(Q5), exact predicates written by us (Q6), own near-body refinement operation
now (Q7), 3D Python prototype with go/no-go (Q8), no penalised CHT (Q9), a
stage 5 for prescribed rigid motion (Q10), signed distance default decided by
benchmark (Q11). The text below has been changed to match.

## 0. Summary

Six stages (0 to 5), each ends in a gallery case with images and a measured accuracy
statement. Every stage keeps the solver data structures untouched: the
immersed body enters as a **mesh operation** (stages 2 to 5) or as a **source
term** (stage 1).

| Stage | Name | What it delivers | Where the wall is | Accuracy class |
|---|---|---|---|---|
| 0 | Geometry kernel | Robust surface checks, winding numbers, signed distance, exact predicates, surface–polygon intersection (D108 part 3) | – | exact geometry |
| 1 | Penalisation pre-stage | Brinkman term on a fixed mesh; density-field topology search; iso-surface to a surface file | smeared over about one cell | laminar and topology only; known flaws listed |
| 2 | Cut cells | Background polyhedral mesh cut by the surface, small cells merged, wall faces are ordinary wall faces, CHT both ways | sharp, on cut cells (wall functions, no layers) | second order laminar; ranking quality for turbulent |
| 3 | Layered distance-field boundary | Thin body-fitted layers grown from a signed distance field, outermost iso-surface cut and imprinted (D58) | sharp, on regular layer cells | target: body-fitted quality at the wall |
| 4 | Hand-off to a body-fitted mesh | The measured immersed-versus-body-fitted check of D211, automatic in the study layer | body-fitted | reference |
| 5 | Prescribed rigid motion | Bodies that translate or rotate through the background mesh: re-cut per step, layers move rigidly, fresh-cell treatment (Q10, 2026-11-22) | as in stage 2 or 3 | as stage 2 or 3 |

| Id | Stage | Content |
|---|---|---|
| D225 | 0 | Geometry kernel: surface check, exact predicates, winding number and signed distance (default by benchmark), tolerant dirty STL, surface–polygon intersection |
| D226 | 1 | Penalisation pre-stage, laminar and low-Re only, density-field topology search |
| D227 | 2 | Cut-cell mesh operation, boundary conditions on cut walls, serial and reproducible |
| D228 | 2 | Small-cell merging, refusal of thin solids |
| D229 | 2 | Wall treatment on cut faces, CHT both ways, own near-body refinement operation |
| D230 | 3 | Layered boundary: heat-type fields, column logic, 3D Python prototype with go/no-go |
| D231 | 3 | Outer iso-surface, imprint join, first-layer loop |
| D232 | 3 | Inward layers for CHT, four-way accuracy comparison |
| D233 | 4 | Automatic hand-off to a body-fitted mesh (the D211 rule as a tool) |
| D234 | 5 | Prescribed rigid motion of cut and layered bodies |

Order is fixed by the maintainer: penalisation first, then cut cells, then
hand-off. Stage 5 (rigid motion) follows the stage it moves (2, later 3).
The layered boundary (stage 3) is the L13 target; it sits between
because it reuses the cutting code of stage 2 and gives the hand-off check
(stage 4) its main use. The hand-off also works after stage 2 alone.

## 1. Existing code this builds on (checked on main at 3359091)

- `foehncraft-mesh::query` (D108 part 1): BVH over surface faces, exact
  closest-point query, deterministic. **Not there:** ray casting, winding
  numbers, signed distance, exact predicates (no orient3d anywhere in the
  crates).
- `imprint.rs` (D58, D182, D186): supermesh by clipping convex polygons,
  closed cells, geometric conservation law. It handles **surfaces of
  revolution only**; arbitrary surfaces need a general surface–polygon
  intersection with exact predicates. This is the shared kernel of stages 2
  and 3.
- Wall treatment D13/D135 and exact wall distance of D48 on any wall face.
- Regions and CHT (D112, D160; 7d1, 7d2 merged): solid cells and interface
  faces exist, so the solid part of a cut cell has somewhere to go.
- Thin walls D189, mesh quality remediation D50, limitations D199, tool
  registry D196, MCP server L18, journal D197.
- **Not there:** the adjoint solver L17 (only D59 kernels), mesh morphing L4,
  adaptive refinement L14.

## 2. Stage 0 – Geometry kernel (D108 part 3, D225)

Scope, in the order of use:

1. **Surface check** tool: closed or open, orientation, self-intersections,
   non-manifold edges, degenerate triangles, duplicate vertices, bounding
   box, units; a report, no silent repair (diagnostic codes `IMM0xx`).
2. **Exact predicates**: adaptive-precision orientation tests (Shewchuk
   1997), **written by us from the paper** (owner decision Q6, strict clean
   room, no dependency); tested against exact rational arithmetic on
   degenerate and near-degenerate inputs.
3. **Winding numbers and inside/outside**: generalised winding number
   (Jacobson et al. 2013) with the fast hierarchical evaluation (Barill et
   al. 2018) on the existing BVH; deterministic. Tolerates dirty STL, which
   exact-Boolean methods do not: the exact mesh-Boolean methods read for
   this plan need manifold, watertight, non-self-intersecting input
   (Cherchi et al. 2022; Trettner et al. 2022).
   **Dirty STL is tolerated, with a report and no silent repair** (Q2).
   Candidates read as abstracts on the publication page of K. Crane
   (2026-11-22; the full papers are read when this stage starts): the heat
   method for generalised signed distance (Feng & Crane 2024), which
   approximates the distance of geometry with holes, noise or
   self-intersections by completing the shape implicitly, with a claimed
   inside/outside more robust than the generalised winding number and any
   discretisation including grids and tet meshes; winding numbers on discrete
   surfaces (Feng, Gillespie & Crane 2023), which complete open curves and
   label regions; the Laplacian for non-manifold triangle meshes (Sharp &
   Crane 2020), which the heat method needs to run on non-manifold STL; and,
   as an optional repair tool only, untangling of self-intersections without
   changing connectivity (Minarcik, Liu, Crane & Li 2026).
4. **Signed distance**: closest-point distance (exists) with the sign from
   the winding number; narrow-band version on the background mesh nodes.
5. **Surface–polygon intersection**: clip a (non-convex, possibly
   non-planar) cell face by the triangle surface with exact predicates and
   indirect (implicit) intersection points (Attene 2020; Cherchi et al.
   2022 describe the technique). Output: polygon pieces with edge links to
   the surface triangle. This is the new core, shared with stage 2 and with
   the general static imprinting of D58.

Verification: unit tests on spheres, cubes and a bunny-class STL; volume and
area of a cut sphere against the analytic value converging at second order;
brute-force closest point and winding number against the fast versions on
random queries; degenerate cases (vertex on face, edge on edge, coplanar)
exact; bitwise equal on 1, 2, 4 threads (D93).

Benchmark against computer graphics (maintainer rule): narrow-band SDF by
(a) closest-point + winding number (ours), (b) fast sweeping on a sparse grid
(Museth 2017 sparse and parallel fast sweeping; Roosing et al. 2019 scan
conversion on the GPU for CFD mesh generation, floating point, needs a
correctly oriented closed surface), (c) the heat method for generalised signed distance (Feng & Crane 2024,
see item 3); the default method is **decided by this benchmark** (Q11). Metrics:
error against the analytic sphere/torus distance, time per million nodes,
behaviour on dirty input (holes, self-intersections, duplicate triangles;
a ready test set is a part of the stage). The zero level set of each field
is extracted with subgrid marching tetrahedra (Baktash, Gillespie & Crane
2025; full paper read before the stage starts) and compared with the input
surface (two-sided Hausdorff distance, manifold and self-intersection
check), which also tests the extraction of stage 1 and step 3e on known
surfaces (owner decision 2026-11-24, graphics-algorithms assessment Q5).

**Status (#310, part 1).** Built in `foehncraft_mesh::surface`: the
surface check (items 1, IMM001 to IMM009, tool `mesh.surface`), exact
and fast generalised winding numbers (item 3; three-term far field, 4e-3 at
β = 2), the signed distance and segment crossings (item 4); the predicates
(item 2) came with #236. **Part 2 (#315):** the benchmark (default:
closest point with the winding-number sign, Q11; fast sweeping wins only on
flipped triangles, the heat method only at holes) and the exact cut of a
mesh as moments (item 5 as volumes, apertures and wall pieces by Gauss–Green
per dimension, 40 to 2500 times more accurate than a linear level-set cut
on a sphere). Open: polygon pieces with implicit points (for a cut that
builds a mesh) and ties of surfaces lying flat in mesh faces.

Sizing: surface check and predicates 2 PRs, winding numbers 1 PR, SDF 1 PR,
surface–polygon intersection 2 to 3 PRs; code about 30 to 60 min per PR,
CI separate. Gallery: none (core), except a **geometry page** with an image
of a cut sphere and its winding-number field (needs an image by the rule).

## 3. Stage 1 – Penalisation pre-stage (D226)

**Role (decided):** topology search and a cheap first look. Not a second
immersed method. Its iso-surface goes into stage 2.

Method: a per-cell solid fraction or density χ ∈ [0,1] (from D108 SDF
volume fractions, or a design variable) and a Brinkman drag term
−(χ/η) u in momentum (Angot, Bruneau & Fabrie 1999), the topology-optimisation
interpolation of Borrvall & Petersson (2003) for the design variable.
Implementation: an implicit diagonal contribution and a source given as a
**field** like D104 (kernels stay pure functions, D59); k, ε/ω, T are held to
zero velocity/solid values in χ ≈ 1 cells. The coupled and the SIMPLEC
solver need the same term inside the Rhie–Chow coefficient (otherwise
pressure decouples in the solid), a check of its own.

Known flaws, stated so that we do not copy them (maintainer rule):
smeared wall of about one cell; wall shear and heat flux not accurate; wall
functions do not fit; accuracy depends on η and the mesh width, with the
penalisation error vanishing with η for fixed mesh but the stiffness
increasing; leakage through thin walls; an industrial aeronautical CFD
code's published penalisation (arXiv 2404.16132, read for this plan, summary
only) reports pressure well predicted, skin friction not very well predicted,
and diverging leading-edge pressure when its wall model was switched on. We
therefore label the stage **laminar and low-Re plus topology** in the
limitations (D199) and refuse a wall function on penalised walls
(diagnostic).

Open point: choice of η (tied to the mesh width, to a time scale, or fixed);
verified by a sweep, not decided in advance.

Verification (D45 style, analytic and manufactured):

| Case | Reference | Check |
|---|---|---|
| Poiseuille/Couette with the wall inside the domain, wall not on a cell face | analytic | order of convergence in h at η ∝ h² and η ∝ h; wall position error |
| MMS of the penalised equation (source in the solid) on 4 mesh types | manufactured | formal order of the discretisation |
| Stokes cylinder/sphere drag in a channel | published Stokes-limit drag with wall correction (Faxén/Haberman type) | drag error vs h |
| Cylinder Re 20 (Schäfer & Turek 1996) | published drag, lift, pressure difference | error vs the same case on our body-fitted mesh |
| Channel with a density field, topology search of a simple duct bend | none (consistency) | iso-surface at 0.5 reproduces the cut surface within one cell |
| Gradient of drag with respect to χ | finite difference | dual numbers (D59) agree to 1e-8 relative |

Gallery 1: cylinder Re 20 penalised against body-fitted (velocity magnitude
and vorticity side by side, with χ); topology search of a duct bend with the
density field and the extracted surface. Recipe: "put a body in the flow
without meshing it". Cases: `immersed/penalised-cylinder.toml`,
`immersed/topology-bend.toml`.

Iso-surface extraction: marching tetrahedra on a decomposed polyhedral cell
(vertex values from cell values by inverse-distance), benchmarked against
dual contouring for sharp features and against subgrid marching tetrahedra
(Baktash, Gillespie & Crane 2025: a manifold, intersection-free surface from
the edge crossings, several patches per tetrahedron, so thin features of a
coarse density field are kept; owner decision 2026-11-24, graphics-algorithms
assessment Q5). Metrics: manifold and self-intersection check, triangle
quality, area and enclosed volume against the field (stage 3 uses the same).

Sizing: kernel and flow hook 1 PR (about 30 to 45 min code), keys and
diagnostics 1 PR, verification table runs about 1 to 2 hours, topology search
and extraction 1 PR. Gradient check on the D59 kernels needs no L17.

**Status (stage 1 complete, #270 and #291).** Built: primitive and `field`
bodies, the implicit sink in the momentum rows and the Rhie–Chow factor, the
topology search `[immersed.topology]` (design variable γ per cell,
Borrvall–Petersson χ(γ) with q = 0.1, hat density filter with its exact
transpose, optimality-criteria updates with a move limit every `every` flow
iterations, the total-pressure power as objective, one-shot after D206) and
the iso-surface (`foehncraft_mesh::isosurface`, marching tetrahedra on any
polyhedral cell, watertight by edge welding, written as `topology.stl`).
Results against the verification table above:

- Poiseuille channel, walls off the cell faces (16 to 256 cells): the
  effective wall stays within half a cell for every K, its sign set by
  where the wall cuts its cell; the mass-flow error of K = c h² falls at
  about first order, under 0.7 % at 256 cells. K = c h leaks on fine
  meshes.
- MMS of the penalised equations (smooth χ): order 2.00 on hexahedra, 2.05
  on polygons. On triangular prisms and the mixed mesh the velocity error
  stalls at about 6 % **without** the sink as well (and with P = 0, and with
  SIMPLEC): a gap of the 2D prism discretisation, not of the penalisation,
  recorded in the verification gaps.
- Stokes drag: a periodic simple cubic array of spheres against the series
  of Hasimoto (1959) and Sangani & Acrivos (1982) replaces the channel with
  a Faxén correction (an exact periodic reference, no wall correction to
  fit). K = 0.01 h²: +6.2 % at 16³ to +2.6 % at 48³; K = 0.03 h² within
  0.2 % from 24³ on (errors of opposite sign cancelling, not a rule).
- Re 20 cylinder: in #270.
- Topology of a duct bend: 28 designs, power 4.69 to 0.37 W, grey fraction
  0.04; the iso-surface at 0.5 is closed and lies within one cell of the
  cut.
- Gradient: the self-adjoint Stokes sensitivity (μ/K) χ′ |u|² V against
  central differences of the converged power, 3.4e-6 relative on 80×40,
  0.5 % at Re 20 (convection makes the problem non-self-adjoint; the exact
  gradient there needs the L17 adjoint). The dual-number check of the row
  stays with L17 A0 (#271).
- Iso-surface benchmark (`tools/bench/isosurface.py`): dual contouring is
  more accurate on sharp edges (box Hausdorff 0.26 against 0.47 cells), MT
  is guaranteed manifold; subgrid MT equals MT for node values linear on
  the tetrahedra, so it pays only with richer node data (stage 3).

Choice of η (open point above, now decided by the sweep): K = μ/η ≈ 0.01 h²
with a time scale free of the flow; the recipe says so.

## 4. Stage 2 – Conservative cut cells (D227–D229)

Method (decided, step 26 §3.2): the background mesh is intersected with the
surface; fluid parts stay ordinary polyhedral cells (D23); cut cells below a
volume fraction merge into a neighbour; the immersed wall is an ordinary wall
face. Literature read or checked for this plan: Ye et al. 1999, Meyer et
al. 2010, Berger & Aftosmis 2012 (from the step 26 assessment, not re-read
here); Crocker, Dubief & Desjardins (arXiv 1411.1004, full text read): second
order for temperature, **first order for momentum statistics**, mixing
instead of merging for cells below volume fraction 0.1, flux continuity at
the interface imposed implicitly, no wall model, 7 to 10 % cost for the CHT
coupling; Ndiaye et al. 2025 (journal article read, inviscid, first order,
three conservative corrections that are either global, or not always
conservative, or locally conservative). Both are Cartesian-grid methods. We
cut an arbitrary polyhedral mesh, and cut cells become ordinary cells, which
keeps conservation exact without a correction step. That claim is checked by
test (below), not assumed.

Steps:

1. **2a Cut.** For every cell touched by the surface, build the fluid
   polyhedron: split faces by the stage 0 intersection, close the cell with
   the surface patch, keep exact topology (watertight cells). A cell split
   into two disjoint fluid pieces becomes two cells. Wall faces may be
   non-planar; the face area vector and centroid use the triangle
   decomposition like all D23 faces (to be confirmed against the geometry
   module in the first PR). Closure check: sum of area vectors of each
   cell below 1e-14 relative.
2. **2b Small cells.** Merge (owner decision Q3) below a volume fraction threshold (a case key,
   default decided by measurement, the literature uses 0.1 to 0.3) into the
   neighbour with the largest shared fluid face; merged cells are polyhedra,
   which D23 allows. Options measured against each other: merging, mixing
   (Crocker et al.), and state redistribution (explicit-scheme method,
   not suited to our implicit steady solvers, recorded only). D50 quality
   remediation applies to what remains.
3. **2c Wall treatment on cut faces.** The wall function needs the true
   first-cell distance: use the exact closest-point distance (D108) per cut
   cell and per merged cell; y+ varies with cut volume fraction, so the
   y+ map and a per-body statistic go to the limitations (D199). No layers
   on cut meshes. Optional thin prism layer extruded from the cut surface
   (step 26 §3.2) is a stage 3 special case and not built here.
4. **2d CHT.** Cut both ways. The solid part of a cut cell becomes a cell of
   the solid region (D112); the cut face is an ordinary region interface
   (D160). Thin solids (thickness below about two cell widths) are not
   resolved: **refuse** with a diagnostic and point to the thin wall of D189
   (owner decision Q4).
5. **2e Boundary conditions.** The D29 set on the immersed surface through
   named surface patches (STL solids or per-triangle tags); moving wall
   velocity of a body at rest or rotating rigidly with a constant angular
   velocity (frame D55). Real body motion is L4.
6. **2f Background refinement.** Owner decision Q7: **our own near-body
   refinement operation**, not a Gmsh size field. It refines cells by
   distance to the surface (and to the layer stack of stage 3) into
   polyhedral cells with hanging-node faces (allowed by D57), steered by a
   size function; deterministic, serial before partitioning. It is the first
   piece of the refinement kernel of L14 (the adaptive, solution-driven
   version comes later on the same kernel). Verification: volume and face
   area conservation, conservative flux across hanging faces on a smooth
   solution (MMS, second order), quality statistics of the refined cells.
7. **2g Reproducibility.** The cut is a serial mesh operation on the global
   mesh, before partitioning; result bitwise independent of threads and
   ranks by construction (D93).

Verification:

| Case | Reference | Expected |
|---|---|---|
| Cut sphere/cylinder volume and area, closure | analytic | volume and area second order, closure 1e-14 |
| Poisson with Dirichlet on a cut sphere (MMS) | manufactured | L2 order about 2; wall gradient order reported (1 to 2) |
| Annulus conduction, Taylor–Couette with the inner cylinder cut | analytic | second order, torque error vs h |
| Inclined Poiseuille channel, walls not on mesh planes | analytic | velocity and wall shear order |
| Concentric-cylinder CHT, conductivity ratio 9 to 900 | analytic | order, interface flux continuity to round-off |
| Stokes sphere, sphere Re 20 and 100 (published drag) | published | drag error vs the body-fitted mesh with the same cell count |
| Cylinder Re 20 and Re 100 (Schäfer & Turek 1996) | published | same as stage 1 case, now sharp |
| Tilted flat plate, turbulent, wall functions | published skin-friction correlation and our own body-fitted run | Cf error vs y+ spread, plus the cut-cell effect on the log-law |
| Conservation: mass, momentum, energy ledger (D119) | – | closed to round-off, same as a body-fitted run |
| Merge threshold sweep (0.05 to 0.4) | – | error and iteration count vs threshold; sets the default |
| Shape-gradient check (§7) | finite difference | measured, including noise from topology switches |

Gallery 2: sphere or cylinder in a duct with the cut mesh closeup and the
velocity field; a heated pin in a channel (CHT, T field); a slanted plate with
the y+ map. Recipes: "cut a surface into a mesh", "choose the merge
threshold". Cases: `immersed/cut-sphere.toml`, `immersed/cut-pin-cht.toml`.
No four-fold or symbol-like arrangements of bodies in images (maintainer
rule); single bodies or irregular patterns.

Sizing: 2a 2 PRs, 2b 1, 2c 1, 2d 1, 2e 1, verification runs 2 to 3 hours in
all; code about 30 to 60 min per PR. Cut robustness on messy surfaces and
the small-cell statistics on complex CAD are **not sized yet**.

**Status (2a, 2b for primitives #313; STL bodies #322; several faces per
cell pair #329).** Built: the cut of
primitive, `field` and STL bodies on any polyhedral mesh with exact
crossings, the closest-point or dual-contouring wall point, face points on
the clipped faces, the validity-driven merge (faces between two cells
joined only where they lie in one plane, otherwise kept apart: two cells
may share several faces, owner decision 2026-11-26) and the per-cell
fallback to fewer points; the Re 20 cylinder ([gallery 313](../gallery/immersed.md#cut-cylinder))
and the Stokes array of STL spheres, second order on hexahedra
([gallery 322](../gallery/immersed.md#stl-cut)). 2c (#357): wall functions
on cut bodies with the exact first-cell distance, sampled one layer out
along the wall normal, y⁺ per body in LIM011, the half-space body `shape =
"plane"`; the 2° plate against White's correlation and the body-fitted plate
([gallery 357](../gallery/immersed.md#cut-plate)). 2d (#359): conjugate
heat transfer cut both ways (the solid part of each cut cell a cell of the
solid region, shared cut faces as the D160 interface), thin solids refused
with IMM205; concentric cylinders second order for k ratios 9 to 900 and a
heated pin in the Schäfer–Turek channel ([gallery 359](../gallery/immersed.md#cut-cht)).
2e (#360): named surface patches of STL bodies with their own conditions,
turning walls with `omega` and their torque, inverted bodies; Taylor–Couette
with both cylinders cut (velocity second order, torque 0.23 % on 176 cells
across) and a cylinder blowing through a slot
([gallery 360](../gallery/immersed.md#cut-patches)). Open: the force balance
on tetrahedra, non-planar polyhedral faces (flat solid slivers on
polyhedra), the k overproduction of the first cut layer, patch edges inside
a wall face, and 2f.

## 5. Stage 3 – Layered distance-field boundary (the L13 target, D230–D232)

As recorded in L13 (PR #173): layers grown from a signed distance field at
program start, the outermost iso-surface as the smooth cut, joined to the
background by imprinting (D58), layers also inward for CHT. The 2D prototype
(step 26 §9) held up except at sharp convex edges and column meeting points;
the heat-field statistics on the fins were unreliable and are not used.

Work packages, each with its own go/no-go (a research stage: results decide
the next one):

1. **3a 3D prototype in Python (numpy/scipy), meshes only; go/no-go before
   any Rust (owner decision Q8).** Same method as
   the 2D prototype, on a sphere, a box, a fin with a narrow gap, an STL with
   a sharp edge. Questions: fan at convex edges (the 2D fan shrank 3.7 times,
   not to zero), merge/stop of columns at concave corners and gaps, local
   layer limit from the gradient dip, cost at scale (the 2D grid solves took
   10 to 15 s per 500k cells). Output: go/no-go for 3D in Rust. About 1 to 2
   hours of agent time per iteration; the number of iterations is not sized.
2. **3b Fields.** Compare for growth directions and layer limit: exact
   closest-point direction (corner kink), heat method (Crane et al. 2013),
   signed heat method (Feng & Crane 2024, approximate, about first order, no
   analysis of corners or medial axis in the paper as read; works on dirty
   geometry, which ties it to the dirty-STL decision Q2). Values for the
   turbulence wall distance stay exact (D108).
3. **3c Surface under the layers.** Surface mesh quality: isotropic remeshing
   (Botsch & Kobbelt 2004), or a quad-dominant parameterisation (Corman &
   Crane 2025, the best candidate named in step 26 §8.6), or the footprint of
   the background mesh on the surface. Measured by the quality of the grown
   layer cells.
4. **3d Column logic in Rust.** Insertion at convex edges, merging and
   stopping at concave corners and gaps, layers down to a minimum (below the
   first-layer height fall back to stage 2 cut cells with a limitation).
5. **3e Outer iso-surface.** Marching cubes against dual contouring and
   marching tetrahedra on the field, and subgrid marching tetrahedra
   (Baktash, Gillespie & Crane 2025: manifold, intersection-free mesh from
   edge crossings, several patches per tet); the field is smooth there, so
   the benchmark measures curvature error, triangle quality and robustness.
6. **3f Join.** Imprint with D58 (general version from stage 0 kernel),
   background size matched to the outer layer height (size jump limit).
   Exactness: conservation as in stage 2.
7. **3g First-layer height loop.** Start from a flat-plate estimate; after
   the first run regrow from measured wall shear (meshing is automatic).
8. **3h Inward layers for CHT.** Negative distances, conformal interface
   with layers on both sides.

Verification: the stage 2 cases repeated on the layered mesh, with the
**four-way comparison** penalised / cut / layered / body-fitted at equal
cell count as the central accuracy statement; layer-mesh quality
metrics (skewness, non-orthogonality, volume ratio); y+ spread; mesh
independence of the immersed surface (the outer cut moved by one cell must
not change the wall result beyond the stated tolerance: this is the
test that the wall really sits on the layers, not on the cut).

Gallery 3: finned heat sink with the layers in the fin gaps (cut view of the
mesh and temperature) — this ties to the heat sink study already running as
its own thread; a sphere with the layer stack and the wake; the four-way
comparison figure. Recipes: "grow layers on an STL", "tune the first layer".

Sizing: 3a not sized yet; 3b to 3h code about 30 to 60 min per PR once 3a
says go; runs 1 to 2 hours; 3D robustness on arbitrary CAD **not sized yet**.

## 6. Stage 4 – Automatic hand-off to a body-fitted mesh (D211, D233)

The rule of D211 as a tool, in the study layer (D62), not in the solver.

1. **`immersed.compare`**: run the same design on the immersed mesh and on a
   body-fitted mesh (Gmsh script of D62 with layers, converted to polyhedra by
   D54), report the difference of each objective, wall quantities and
   y+, with the grid-convergence verdict (D67). Record as a limitation
   (D199). The difference becomes the first-order correction of D207 (the
   immersed mesh is a fidelity level, tagged in D106).
2. **Switch rule**: switch when the trust-region radius or the gain per round
   falls below the measured difference, or the gain per round falls below the
   uncertainty of the immersed result (D206). "No switch needed" is a valid
   outcome, especially after stage 3.
3. **Surface export** for the body-fitted mesh: final design surface to STL,
   plus a template with the meshing parameters. Morphing (L4) instead of
   regenerating is later.
4. **Agent safety**: starts runs, so it is a `start-run` tool under the D198
   budgets; an estimate comes first (D81).

Verification: the comparison is itself the test; it must reproduce, on the
sphere, cylinder and finned plate, the differences that stage 2 and 3 report.

Gallery 4: a small fin-pitch study on the immersed mesh, then the check on the
body-fitted mesh; plot of objective versus round and the temperature fields of
both meshes. Recipe: "study on cut cells, check on a fitted mesh".

Sizing: tool and report 1 PR (about 30 to 45 min code); study runs 1 to 2
hours.

## 6a. Stage 5 – Prescribed rigid motion (owner decision Q10, D234)

Added on 2026-11-22. Bodies translate or rotate along a prescribed path
(expression D53 or table; constant angular velocity as a special case)
through the background mesh. Deformation and free (6-DOF) motion stay with
L4; this stage only moves rigid bodies.

- **Cut cells (stage 2):** the cut is redone every step, or when the body has
  moved by a set fraction of a cell (a case key). Cells that change from
  solid to fluid (**fresh cells**) get velocity, pressure and scalars from
  their neighbours by a defined extrapolation, and the step reports how many
  cells changed. The literature on the fresh-cell problem is read when the
  stage starts (not read for this plan).
- **Layered boundary (stage 3):** the layer mesh moves rigidly with the body
  and is not rebuilt; only the outer cut and its imprint (D58) change per
  step. This reuses the imprinting of the sliding interfaces (D182, D186)
  generalised in stage 0.
- **Conservation:** the geometric conservation law must hold for the faces
  that appear, vanish and move (as in D186); mass and energy ledgers (D119)
  close to round-off.
- **Merging** of small cells is recomputed per step; the merge map must not
  flip between steps in a way that makes the solution jump (check by a
  smooth-motion test).

Verification:

| Case | Reference | Check |
|---|---|---|
| Cylinder dragged through still fluid against the same cylinder fixed in uniform flow | Galilean equivalence, our own fixed-body run | drag and lift history agree within the stage 2 error |
| Sphere in accelerated motion in a still fluid, low Re | analytic added-mass and history force limits (to be sourced when the stage starts) | force history |
| Rotating body: Taylor–Couette spin-up with the inner cylinder cut and rotating | analytic, and our sliding mesh (D182) | torque history |
| Smooth sub-cell motion of a body | none (consistency) | no jump in forces when a cell changes state; merge-map flips counted |
| Ledgers (D119) | – | closed to round-off including fresh cells |

Gallery 5: a body of 3 or 5 blades (no four-fold pattern) turning through the
background mesh, velocity and vorticity; a cylinder dragged through still
fluid with its wake. Recipe: "move a body through a mesh without remeshing".

Sizing: cut per step 1 PR, fresh-cell treatment 1 to 2 PRs, conservation
checks 1 PR, code about 30 to 60 min per PR; the whole stage **not sized
yet** because the fresh-cell treatment is research.

## 7. Interaction with the adjoint (L17) and CHT

**The adjoint solver is not built.** All gradient statements below are plans
for when L17 exists; stage 1 to 4 do not wait for it.

| Method | What the adjoint gives | Shape gradient | Risk |
|---|---|---|---|
| Penalisation | exact discrete gradient with respect to χ per cell (D59 kernels already agree with finite differences) | topology gradient directly | none for gradients; flaws are in the wall accuracy |
| Cut cells | discrete adjoint on the frozen cut mesh | (a) continuous surface sensitivity on cut wall faces (the L17 surface sensitivity), applied as normal displacement; (b) forward-mode derivative of the cut geometry itself with frozen topology, generic over the D7 float type so dual numbers pass through the cutting code | cut topology and merge decisions change discontinuously with the surface; the gradient is piecewise smooth, noisy at switches |
| Layered | layer nodes are smooth functions of the surface and design parameters | adjoint node sensitivities of L17 chained with the layer-mesh Jacobian; outer cut far from the wall treated as frozen with an error estimate for D206 | needs the column logic to be differentiable (frozen topology, smooth rules) |

Rule for stage 2: the gradient strategy is decided by measurement when L17
exists, not in advance: compare (a), (b) and finite differences (a step sweep
shows the noise) on the cylinder and the pin. The cut code is written
generic over the float type from the first PR, so (b) stays possible. A
shape-and-topology derivative on a fixed mesh with a level set exists in the
literature (Gangl & Gfrerer, arXiv 2209.15491, read as a summary: 2D, linear
level set, mesh fixed) and informs the pre-stage; it is not adopted for the
3D design chain.

**CHT.** Cut and layered methods cut both ways; the interface is a normal
D160 face with series conductance, so interface fluxes are conservative and
the thermal adjoint (D160 kernels) applies unchanged. Penalisation: solid
temperature through an equation domain of regions (D112) is not available,
because penalised CHT is left out of stage 1 (owner decision Q9). Cycle-averaged CHT (D204) works on top of either mesh.

## 8. Interfaces: scriptable and MCP compatible (D77–D87, L18, D196)

Nothing exists only in a GUI. New tools in the D196 registry, typed in/out,
action kind per D198:

| Tool id | Kind | Does |
|---|---|---|
| `geometry.check` | read | surface report (closed, orientation, self-intersection, volume, bounds) |
| `geometry.sdf` | write-files | signed distance and winding number on a mesh or grid, VTKHDF |
| `mesh.cut` | write-files | background mesh + surface → cut mesh, quality and small-cell report |
| `mesh.layers` | write-files | surface → layered mesh with the outer cut surface |
| `immersed.estimate` | read | cell counts, small-cell count, layer count, memory before the run (D81) |
| `immersed.compare` | start-run | the D211 check |
| `study.immersed` | start-run | the study with automatic hand-off |

Case file keys (D77 schema, Python mirror D61):
`[immersed] method = "penalisation" | "cut" | "layered"`, `[[immersed.body]]`
with `surface`, `units`, `kind = "solid" | "region"`, `material`, patch tags;
`[immersed.cut] min_volume_fraction`, `merge = "neighbour"`;
`[immersed.layers] first_height | y_plus_target`, `growth`, `count |
total_thickness`, `min_layers`; `[immersed.penalisation] eta`, `interpolation`.
Diagnostics `IMM001`… with fix hints (D78): surface not closed, self
intersection, body outside the domain, thin solid, first-layer height below
the cut, small-cell fraction too high. Limitations (D199): number and
fraction of merged cells, y+ spread on immersed walls, thin-wall fallbacks,
penalised walls. Machine-readable report in `result.json`. Every tool is
deterministic and journaled (D197). Surface inputs: STL (ASCII, binary),
OBJ later; STEP through Gmsh/OpenCASCADE (D62). Verification grades go into
`features.toml` per stage (D199): stage 1 "partial, laminar", stage 2
"verified" only after the table above passes.

## 9. Owner decisions (2026-11-22)

All eleven questions were answered by the maintainer on 2026-11-22, one card
at a time. Each line gives the question, the answer and what it changed.

| # | Question | Answer | Effect on the plan |
|---|---|---|---|
| Q1 | Accuracy targets (proposed 1 % laminar drag, 5 % turbulent wall shear, 2 % for stage 3) | **Set no numbers yet** | Stage 2 reports measured errors only; targets are fixed after the first four-way comparison of stage 3 |
| Q2 | Dirty STL: tolerant or closed manifold only | **Tolerant, with a report** (maintainer: use Crane's methods, for example the latest heat method, for dirty geometry) | §2 item 3: the heat method for generalised signed distance, winding numbers on discrete surfaces, the non-manifold Laplacian and optional untangling become candidates; no silent repair |
| Q3 | Small cell handling | **Merge**; mixing and state redistribution measured as references only | §4 step 2b |
| Q4 | Solids thinner than about two cells | **Refuse**, point to the thin wall of D189 | §4 step 2d |
| Q5 | Penalised walls with turbulence | **Laminar and low-Re only**, wall function refused | §3 |
| Q6 | Exact predicates | **Write them ourselves** from the paper | §2 item 2 |
| Q7 | Background refinement near the body | **Own refinement operation now** (instead of a Gmsh size field) | §4 step 2f, first piece of the L14 refinement kernel |
| Q8 | 3D Python prototype of the layered boundary first | **Yes**, go/no-go before Rust | §5 step 3a |
| Q9 | Penalised CHT in stage 1 | **Left out**; CHT starts at stage 2 | §3, §7 |
| Q10 | Moving bodies | **Plan a stage for prescribed rigid motion** | new stage 5, §6a |
| Q11 | Default signed distance method | **Decided by the benchmark** of §2 | §2 benchmark |

The ids D225 to D234 are entered in requirements.md with this record.

## 10. Order of work and what stops what

Stage 0 and stage 1 can start together (no overlap in code). Stage 2 waits
for stage 0 steps 1 to 5. Stage 3a (Python) can run next to stage 2 and
needs only stage 0 steps 3 and 4. Stage 4 waits for stage 2 and uses stage 3
when it exists. Stage 5 waits for stage 2 (and for stage 3 to move layers);
the near-body refinement operation of step 2f is built with stage 2 and is
useful on its own. Nothing here depends on L17; the shape-gradient checks of §7
wait for it. Every PR: changelog fragment, README "What works today", docs,
gallery entry with at least one image and field images for physics (the
maintainer rule), case files and a recipe.

## 11. Sources

Read for this plan (full text or the part named): Crocker, Dubief &
Desjardins, arXiv 1411.1004 (full); Ndiaye et al., Int. J. Heat Fluid Flow
114, 2025 (journal article, part); Cherchi, Pellacini, Attene & Livesu,
arXiv 2205.14151 (part); Trettner, Nehring-Wirxel & Kobbelt, ACM TOG 41(4),
2022 (part); Roosing, Strickson & Nikiforakis, arXiv 1903.00353 (part);
arXiv 2404.16132 (part, summary); Gangl & Gfrerer, arXiv 2209.15491 (part,
summary); Bayat, Li & Alexandersen, arXiv 2601.02202 (abstract level:
density-based topology optimisation for turbulent flow with a Brinkman
penalty and an implicit wall function).
From the earlier step 26 assessment (not re-read now): Ye et al. 1999;
Meyer et al. 2010; Berger & Aftosmis 2012; Peskin 2002; Mittal & Iaccarino
2005; Angot et al. 1999; Borrvall & Petersson 2003; Tseng & Ferziger 2003.
Read as abstracts on K. Crane's publication page (2026-11-22): Feng &
Crane 2024 (heat method for generalised signed distance), Feng, Gillespie &
Crane 2023 (winding numbers on discrete surfaces), Sharp & Crane 2020
(Laplacian for non-manifold meshes), Baktash, Gillespie & Crane 2025
(subgrid marching tetrahedra), Minarcik, Liu, Crane & Li 2026 (untangling
surfaces). Named from memory and to be read when their stage starts: Shewchuk 1997;
Jacobson et al. 2013; Barill et al. 2018; Attene 2020; Museth (sparse fast
sweeping); Botsch & Kobbelt 2004; Crane, Weischedel & Wardetzky 2013;
Schäfer & Turek 1996. A "continuous adjoint cut-cell method for cavitating
flows" turned up in the search by title only and is unread.
