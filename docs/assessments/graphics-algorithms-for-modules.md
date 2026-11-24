# Computer-graphics and geometry-processing algorithms for the Foehncraft modules

Assessment, 2026-11-24. Requested by the maintainer: "I need the best modules
and want to speed them up with the computer graphics algorithms, for example
Keenan Crane's radiation and Monte Carlo algorithms." Discussion only, no
code. It applies the standing rule of 2026-11-21 (for geometry and interface
topics always benchmark computer-graphics algorithms; the goal is
state-of-the-art accuracy).

It does **not** repeat what the earlier pages already decided. Read first:
the [spray prerequisites](spray-module-prerequisites.md), the
[L13 plan](../requirements/l13-immersed-boundary-plan.md), the
[MPI sizing](mpi-multi-rank-sizing.md), and in the project files
`radiation-monte-carlo-assessment.md` (L15; D107, D108),
`design-optimization-surrogates-immersed-assessment.md` (section 8, layered
distance fields) and `wall-film-assessment.md`. Where this page adds to one
of them it says so.

**What was read for this page.** The code on `main` at 84bf123 (the query,
locate, renumber and partition modules, the exact sum), the requirements
rows L12, L13, L15, D48, D107, D108, P9, the six assessments above, the
current publication list and abstracts on Keenan Crane's page
(cs.cmu.edu/~kmcrane, read 2026-11-24, because earlier work was found to be
newer than memory), and the licence files of the libraries named in
section 8. **Not read in full:** the papers themselves. For every Crane
paper only the abstract on his page was read, and for the other papers only
the bibliographic data. That is enough to rank them; the full text is read
when the stage that uses it starts (as the L13 plan already says). No source
code of any library was read, only licence texts. Items found by search but
not read are listed in section 9.

## 1. Answer in short

1. **The largest sure gain is small and unglamorous: a hierarchy over cell
   boxes and exact geometric predicates in the geometry-query layer** (D108
   part 2 and the predicates of L13 stage 0). They serve point location,
   mesh-to-mesh interpolation, probes, particle tracking at the rare places
   where the walk fails, and ray casting. Everything else stands on them.
2. **Radiation is where graphics algorithms give the largest speed and
   accuracy gain, but not through the ray-tracing hardware or a faster
   BVH.** In a gas mesh the bundles walk from cell to cell and need no
   hierarchy. The gain comes from volume path-tracing methods: null-collision
   tracking on a coarse majorant grid (cost proportional to the number of
   collisions, not to the number of cells crossed), stratified and
   low-discrepancy emission, ratio tracking for the next-event estimates and
   multiple importance sampling between the forward and reciprocal
   estimators. The first three keep our bitwise reproducibility; the
   majorant grid needs care at rank boundaries.
3. **Crane's heat method does not improve the RANS wall distance.** The
   exact distance by BVH closest point (already on main, D48) is exact at
   the first cell, deterministic and differentiable almost everywhere; the
   heat method is approximate, with an error that depends on mesh size and
   time step. The heat method family earns
   its place where the geometry is *dirty or implicit*: signed distance and
   inside/outside for STL bodies (L13 stages 0 and 3), which the L13 plan
   already benchmarks.
4. **Walk on spheres and its relatives stay what the radiation assessment
   said:** a verification tool and a thick-region step inside the photon
   Monte Carlo, not a production solver for the flow, the film or the solid
   (section 5).
5. **Sprays, wall films and the solver's mesh operations gain little from
   fancy graphics methods.** The walk already is the right tracker; the
   film is a finite-volume problem on a polyhedral surface; the real
   gains there are exactness (predicates) and cheap ordering and sampling
   tricks (section 6).
6. **Do not use the ray-triangle test known as "watertight"** (Woop, Benthin
   and Wald 2013) without a patent check: patents with that title exist
   (section 8). A watertight test built on our own exact orientation
   predicates is simple and avoids the question.

Recommended order: section 10. Questions for the maintainer: section 11.

## 2. State of the code that these methods would change

| Item | State on main | Note for this assessment |
|---|---|---|
| BVH over wall faces, closest point | `query/bvh.rs`: binned SAH build (16 bins, leaves of a few triangles), depth-first query, nearer child first, ties by smallest face id | Exact wall distance D48 is done; the build follows MacDonald and Booth 1990 and the binned form of Wald 2007 |
| Ray casting, cell-box hierarchy, winding numbers | not there (D108 parts 2 and 3) | Needed by radiation surface exchange, film impact on thin walls, mesh-to-mesh interpolation, L13 |
| Exact predicates | not there (`foehncraft-core::exact` is the exact *sum* of D93, not geometry) | L13 stage 0 / D225 writes them |
| Point location | `locate.rs`: tetrahedral-decomposition walk (Macpherson, Nordin, Weller 2009) with plain floating-point signed volumes; if the walk is blocked, a **uniform grid of cell boxes** with bin size twice the mean cell size, at most 1024 bins per axis, and a linear test of the cells in the bin | The fallback is the weak point (section 3.1). It is also the whole algorithm for any query without a good hint cell |
| Cell order | `renumber.rs`: reverse Cuthill-McKee; random order for the worst case | G1 measured 3 to 5 times less bandwidth for a random order |
| Partition | METIS (feature), recursive coordinate bisection | |
| Reproducibility | exact long-accumulator sums, counter-based RNG with stream words, fixed tie-breaking (D93) | Every method below is judged against this |
| Float width | float64 everywhere, P9 (single-precision traversal and vendor ray-tracing units) is an option outside the strict mode | |

## 3. The geometry-query layer (D108 parts 2 and 3)

### 3.1 Hierarchy over cell boxes replaces the uniform bin fallback

**Observation (from reading `locate.rs`, not measured).** The bins are sized
by the *mean* cell size. On a graded mesh (boundary layers, near-wall tets,
refinement zones) a single bin then holds many small cells, tested one by
one, and a large cell is listed in many bins. The walk from a good hint cell
is the primary path, so particle tracking is hardly affected; but **mesh-to-mesh
interpolation (D52), probes and line samples, restart on a changed mesh and
the first location of injected parcels have no good hint** and go straight
to this search.

**Method.** A BVH over the cell bounding boxes (the same binned-SAH build as
on the wall faces, so no new code path), query "which cell contains x" by
descending into every box that contains x and testing the cell with the
tetrahedral decomposition of `locate.rs`. Standard graphics practice; the
survey is Meister et al. (2021). The result is O(log n) independent of the
size grading.

| | |
|---|---|
| Module | point location, D52 interpolation, probes, injectors, AMR donor search (L14) |
| Expected gain | removes a size-grading dependence; speed on graded meshes to be measured (benchmark: all cell centres of the heat exchanger mesh as queries, grid fallback against BVH, with and without a hint) |
| Cost | one PR, about an hour of code plus one CI cycle |
| Risk | low; build is deterministic if the sort uses total orders only (as `bvh.rs` already does) |
| Bitwise reproducibility | kept: when two cells claim a point on a shared face or edge, the smallest global cell id wins, as the wall query does |
| MPI | the hierarchy is built per rank over its cells and halo; donor search across ranks needs the global box hierarchy of the rank boxes (small), the same design as the replicated wall BVH |

### 3.2 Exact geometric predicates and a tracker that cannot lose a parcel

The tracker decides the exit face from signed tetrahedron volumes in
floating point. Where a track passes exactly through a vertex or an edge
(common on structured and extruded meshes) the signs can disagree between
the two cells of a face by rounding, and a parcel is lost or loops. The
spray assessment already lists "lose no parcel silently" as a risk.

**Method.** Shewchuk's adaptive-precision orientation test (Shewchuk 1997):
a floating-point evaluation with a static error bound, an exact fallback
only when the sign is uncertain. The L13 plan (D225, Q6) already has us
write these predicates from the paper. **The addition of this page:** use
them in the tracker too, as a filtered test, so the exit decision is exact
and identical in both cells of a face. Where a *constructed* point is
needed (the crossing point itself), indirect predicates (Attene 2020) avoid
rounding the construction.

| | |
|---|---|
| Module | spray tracker, radiation bundle tracker (same walk), L13 stage 0 and 2 |
| Expected gain | correctness, not speed: zero lost parcels and bundles at degenerate crossings, and robustness for the dirty-STL work |
| Cost | none beyond D225; the filter makes the common case about as fast as the plain test (to be measured on the tracker tests) |
| Risk | low; the predicate is deterministic by construction, so it is *more* reproducible than rounding-dependent signs |
| Reproducibility, MPI | kept and improved: both ranks that share a processor face decide the same way |

### 3.3 Ray casting on the BVH, and what to build it from

Needed by: surface-to-surface radiation and view factors, spray impact on
walls that are not mesh faces (thin walls of D189, immersed bodies), film
impact points, headless renders (L21), and the inside/outside ray tests.

**Method choices, in order of preference.**

1. **Ray-triangle test from our own exact orientation predicates** (the
   sign of the three edge orientations decides a hit; consistent on shared
   edges, so the test is watertight by construction and uses nothing but
   the predicate of section 3.2). Möller and Trumbore (1997) is faster but
   not watertight; use it only as the non-strict fast path with an exact
   fallback when a sign is close to zero.
2. **Robust traversal** of the box test: conservative rounding of the slab
   test (Ize 2013) so that a ray never slips between two boxes. In float64
   with exact box bounds the issue is much smaller than in the float32
   renderers the paper addresses; keep the guard anyway because it is a
   one-line multiplication.
3. **No epsilon offsets for secondary rays.** Renderers use offsets against
   self-intersection (Wächter and Binder 2019 give an error-bound form).
   Our rays start at known mesh points, so start a ray *inside a known cell*
   and let the cell walk carry it; the hierarchy is then only needed for
   surfaces that are not mesh faces. This removes the problem.

| | |
|---|---|
| Expected gain | enabling, not a speed-up |
| Cost | two PRs (closest hit and any hit, then the batched GPU form), not sized beyond that |
| Risk | the watertight test of Woop, Benthin and Wald (2013) is widely used but patents with the title "Watertight ray triangle intersection" exist (section 8); our own predicate-based test avoids it |

### 3.4 Faster hierarchies: wide nodes, SIMD, GPU build

The state of the art is well documented (survey: Meister et al. 2021):

- **Wide BVH (4 or 8 children per node) with SIMD box tests** (Dammertz,
  Hanika and Keller 2008); on the GPU, compressed wide nodes (Ylitie, Karras
  and Laine 2017). Published gains are a small multiple of the binary tree
  on float32 rendering workloads; for float64 CPU the SIMD width is half and
  the gain is smaller. **Not measured for us.**
- **GPU build:** linear BVH from Morton codes (Karras 2012, already named in
  D108), improved by treelet restructuring (Karras and Aila 2013) or parallel
  locally-ordered clustering (Meister and Bittner 2018). Refit for moving
  walls is the cheap operation and is already in D108.
- **Queries for all cell centres at once** (the wall distance sweep): sort
  the queries along a Morton curve so neighbouring queries share traversal
  paths, and start each query with the bound `d(x) <= d(y) + |x - y|` from
  an already finished neighbour `y` (triangle inequality), so whole subtrees
  are skipped from the first node. Standard in distance-field computation;
  no new data structure.

| | |
|---|---|
| Module | wall distance, ray casting, radiation traversal where used |
| Expected gain | none needed for the wall distance today (the radiation assessment estimated seconds on a GPU for 100 M cells, one pass, rerun only when walls move); worth doing only when the per-call ray counts of radiation or the moving-wall refit show a need |
| Cost | wide BVH one PR; GPU build one to two PRs, not sized |
| Risk | float64 SIMD gain may be small; the GPU build needs the D60 backend kernels |
| Reproducibility | the tree shape may differ between builds if a parallel build is non-deterministic: use total orders and fixed tie-breaks (as the host build does); the query result is independent of tree shape as long as ties go to the smallest id |

**Recommendation:** measure first. The present wall distance is a one-off
per mesh; nothing in the current modules is limited by it. Do 3.4 when the
radiation ray counts exist.

## 4. Radiation (L15): what the rendering literature adds

The radiation assessment fixed the design (photon Monte Carlo for the RTE,
reciprocal estimators, counter-based RNG, fixed-point tallies, sphere steps
in thick regions). The additions below are the rendering methods *not yet
in it*, ranked by value over cost.

### 4.1 Where the time goes, and the one structural choice

A bundle in a gas mesh advances cell to cell with the tetrahedral walk
(D47). The cost per bundle is the number of cells crossed, which in a
thin gas (long mean free path, fine mesh) is hundreds. **Null-collision
(delta, Woodcock) tracking on a coarse majorant grid** replaces this: a
coarse grid stores the maximum absorption coefficient per super-cell; a
free flight is sampled against that majorant, and the *real* coefficient is
looked up only at the sampled point (by point location, section 3.1);
a collision is real with probability `k(x) / k_majorant`. The cost is the
number of collisions plus the number of super-cell crossings, not the number
of cells (Woodcock et al. 1965 for the method; Galtier et al. 2013 and
Eymet et al. 2013 for its use in radiative transfer;
Szirmay-Kalos et al. 2011 and Kutz et al. 2017 for the rendering form with
hierarchical majorants).

| | Walk through cells | Null-collision on majorant grid |
|---|---|---|
| Cost per bundle | cells crossed | collisions + super-cell crossings |
| Thin gas, fine mesh | high | low (hypothesis: up to an order of magnitude; to be measured) |
| Thick, strongly varying gas | cheap (few cells) | many null collisions where the majorant is loose |
| Absorption tally | exact path length per cell (track-length estimator) | collision estimator (noisier) unless the walk is also done |
| Ledger closure (D107) | exact per bundle | exact per bundle (every bundle has a fate) |
| Rank boundary | stops at the processor face | super-cells must stop at the rank boundary, or the bundle hands over mid-flight |
| Risk | none | medium: noisier tally, handover logic |

**Recommendation:** build R1 first on the walk (it is the reference and
the exact tally), then add the majorant-grid tracker as a switchable
estimator and **benchmark both on the heat exchanger and a thin-gas furnace
case** (accuracy per cost: same noise level, time per radiation call). The
decision between them is by that benchmark, as the maintainer's rule asks.
Both keep the bitwise behaviour: the random draws are keyed by (emitter,
bundle, call, event) in either; the number of draws per bundle is a fixed
function of the bundle's own history.

### 4.2 Sampling

| Method | What it does | Gain | Cost, risk, reproducibility |
|---|---|---|---|
| **Stratified inverse-CDF emission** | Bundle `i` of `N` takes the uniform number `(i + xi_i) / N` through the cumulative emitted power over the emitters. The selection is exactly proportional to emitted power, has no sampling variance in the *number* of bundles per emitter, and the bundles come sorted by emitter (coherent memory access, a free sort) | removes the emitter-count noise; coherence | trivial; the keyed RNG still gives `xi_i`; the cumulative sum is an exact sum (D93); alias tables (Walker 1977, Vose 1991) are O(1) but lose the stratification and the order |
| **Owen-scrambled Sobol points** for direction, wavelength and position (Burley 2020) | already in D107 / S6 | converges faster than 1/sqrt(N) for smooth integrands | planned; keyed by bundle identity so rank-independent |
| **Ratio and residual-ratio tracking** (Novak et al. 2014) | unbiased transmittance estimates with weights instead of analogue absorption, for the next-event estimate towards walls and for the reciprocal estimator | removes the zero-or-one noise of the transmittance | low; part of the estimator design, same RNG keys |
| **Multiple importance sampling between the forward and the reciprocal estimator** (Veach and Guibas 1995) | the assessment's reciprocal estimator is excellent near isothermal conditions, the forward one far from it; MIS combines them with provably bounded variance instead of choosing a regime | a single estimator that is not bad anywhere | medium: new combination for the net exchange; **a spike of its own**; keep the ledger by booking the combined weight |
| **Hero-wavelength style sampling for the spectrum** (Wilkie et al. 2014) | sample one spectral point (k-value) per bundle and weight the others with MIS | cheaper than looping over gray gases; fits WSGG and FSK by sampling | later, with the spectral stage |
| **Path guiding** (learning where important bundles come from; Müller, Gross and Novak 2017) | radiation is called again every few iterations with slowly changing fields, so a learned importance distribution can be reused | possibly large for noisy cells in cold regions | research; must be deterministic (learn from the fixed-point tallies of the previous call) and unbiased (MIS); do not start before the baseline exists |

### 4.3 GPU form and coherence

Wavefront tracing and bundle sorting by cell (Laine, Karras and Aila 2013) are in the radiation assessment. The rendering literature adds **reordering of rays by the cell or the direction they will visit**, to raise memory coherence (Eisenacher et al. 2013, sorted deferred shading; the idea is the same as sorting spray parcels by cell). With stratified emission the first sort is free. Hardware ray-tracing units (single precision, vendor specific) stay P9. Embree has a permissive licence and a SYCL GPU path (section 8) but is float32; use it as a **timing reference in the benchmark harness**, not as a dependency.

### 4.4 Radiation and the other thermal paths

A published line of work couples conduction, convection and radiation in one
Monte Carlo algorithm on complex geometry (found by search, section 9, not
read). It would matter only if the solid (CHT) were to be solved grid-free
too; our solid is a mesh solve on AMG, which is cheap. Keep as a note.

## 5. Crane's Monte Carlo family and the heat method: where each one fits

| Method | Best use here | Not for | Verdict |
|---|---|---|---|
| Walk on spheres, variable coefficients, walk on stars, Robin boundaries (Sawhney and Crane 2020; Sawhney et al. 2022, 2023; Miller et al. 2024); Zombie (MIT) | R2 sphere steps in thick regions; **verification probes with error bars** for Laplace, Poisson and screened Poisson on complex geometry (S4); thin-feature fields where meshing is the problem | production P1; flow; any non-linear or turbulent problem | as the radiation assessment concluded; unchanged |
| Boundary value caching, mean value caching (Miller et al. 2023; Bakbouk et al. 2023) | reduce variance of many nearby WoS estimates; a candidate for the verification probes | | use if S4 gets built |
| Participating-media WoS (Miller et al. 2025) | not a radiation solver (it estimates the average over random particle configurations) | | parked, as before |
| Differential WoS (listed on Crane's page, abstract not read) | shape derivatives of Monte Carlo fields | | note for L17 only |
| Heat method for distance (Crane, Weischedel and Wardetzky 2013, CACM 2017) | smooth distance-like fields on meshes by two sparse solves (a heat step and a Poisson solve); runs on our AMG | the exact wall distance (approximate; no exactness at the first cell; the error depends on the mesh size and on the time step `t ~ h^2`) | do not replace D48 |
| Vector heat method (Sharp, Soliman and Crane 2019) | parallel transport of vectors and extension of velocities over a surface | | no module needs it now; keep for level-set extension (L3) and surface mesh work |
| Heat method for generalised signed distance (Feng and Crane 2024) | inside/outside and distance for STL with holes, noise or self-intersections | clean CAD (use closest point plus winding number) | already a Q11 benchmark candidate in the L13 plan |
| Winding numbers on discrete surfaces (Feng, Gillespie and Crane 2023) | label regions of open or dirty surfaces | | in the L13 plan |
| Laplacian for non-manifold meshes (Sharp and Crane 2020) | lets the heat method run on raw STL | | in the L13 plan |
| **Subgrid marching tetrahedra** (Baktash, Gillespie and Crane 2025/26) | iso-surface extraction on a tetrahedral background mesh with several surface patches in one cell (thin sheets); also converts polygon soup into a manifold, intersection-free mesh | | **new for us**: candidate for the iso-surface of the penalisation stage (D226), for the dirty-STL tolerance of D225 and for flame-front or interface surfaces in post-processing; add to the L13 benchmark list; abstract only read |
| Intrinsic triangulations, integer coordinates, intrinsic simplification (Sharp, Soliman and Crane 2019; Gillespie, Sharp and Crane 2021; Chislett et al. 2023) | well-conditioned surface Laplacians on bad triangulations; simplify an STL while keeping the Laplacian; geometry-central is MIT | | useful *inside* the signed heat method on rough STL; not before |
| Rectangular surface parameterization (Corman and Crane 2025) | quad-dominant surface mesh under the layers of L13 stage 3 | | already noted in step 26 |
| Points as tori, ray tracing harmonic functions (Feng, Gkioulekas and Crane 2026; Gillespie et al. 2024) | signed distance for point clouds; sphere tracing of harmonic functions | | no module; the screened-Poisson theory is background reading |
| Monte Carlo fluids (Rioux-Lavoie et al. 2022; velocity-based follow-up by Sugimoto, Batty and Hachisuka) | grid-free vorticity or velocity methods for graphics-quality fluids | RANS, wall functions, energy, compressibility, industrial accuracy | **not recommended**: no turbulence closure, no wall model, noise in every cell; abstracts and titles only |

**Why the heat method does not help the RANS wall distance.** D48 needs
the exact Euclidean distance to the nearest wall point at the cell centre,
exact at the first layer, because the near-wall models are calibrated with
it (y-plus at the first cell). The BVH answer has all of that, returns the
wall face and normal, and its derivative is the unit vector to the closest
point (almost everywhere), which is all L17 needs. The heat method gives a
smooth field with a discretisation error that grows with the mesh size and
the time step; the L13 plan already treats it as one of three signed-distance
candidates, which is the right place for it.

## 6. The other modules

### 6.1 Sprays (L12)

| Idea from graphics | Use | Verdict |
|---|---|---|
| Cell-box BVH, exact predicates (sections 3.1, 3.2) | first location of injected parcels, lost parcels at degenerate crossings | do |
| **Stratified (quantile) sampling of the droplet size and the cone** instead of independent random draws | a Rosin-Rammler, log-normal or table distribution is represented by `N` parcels at the quantile midpoints (or stratified with the keyed random number inside each stratum), which should cut the parcel count needed for a given fidelity of the size distribution (expected, to be measured); Sobol points for the cone directions | do, with SP1; keep the keyed RNG so the choice is rank independent |
| Sort parcels by cell (Morton or cell id) every step | coherent memory access and the sorted deposit that D117 already requires | already implied; use the cell order of `renumber.rs` |
| Ray casting for wall impact | only for walls that are not mesh faces (thin walls, immersed bodies) | later, section 3.3 |
| Spatial hashing, sweep and prune for collisions | the O'Rourke and Nordin collision models are per-cell stochastic and need no pair search | not needed |
| Space-filling-curve partition for the spray load | repartition by cutting a sorted list (cheap, incremental) | the halo and METIS decide the mesh partition; keep as an option for SP3 only |

### 6.2 Wall film (L12 film, D123, D124)

The film lives on the surface mesh of the wall faces (polyhedral, finite
volume). Parallel transport of the film velocity between neighbouring faces
is a rotation about the shared edge by the dihedral angle: exact and trivial
on a polyhedral surface, so the vector heat method is not needed. Intrinsic
triangulations would matter only if the film moved to a finite-element
Laplacian on a poor triangulation. **Verdict: nothing from graphics is
worth building for the film now.** The impact point comes from the tracker
(exit face), not from a ray cast. Revisit if a film on immersed (STL)
surfaces is wanted.

### 6.3 Mesh and interface operations

| Operation | Graphics method | Verdict |
|---|---|---|
| Sliding and non-conformal interfaces (D58, D182) beyond surfaces of revolution; cut cells; imprint of STL | exact predicates and indirect constructions (Shewchuk 1997; Attene 2020); mesh arrangements for robust Booleans (Cherchi et al. 2022, which need clean input) | already the L13 plan; use the same kernel for general D58 imprint |
| Candidate pairs for supermesh clipping | BVH over face boxes (section 3.1 code) | do with D108 part 2 |
| Cell ordering for caches | Hilbert or Morton order (cache-oblivious layouts: Yoon et al. 2005) against reverse Cuthill-McKee | a small benchmark in the G1 style; expected gain small; RCM stays unless the benchmark says otherwise |
| Mesh-to-mesh interpolation (D52, restarts, AMR L14) | cell-box BVH for the donor, walk to finish | section 3.1 |
| Surface repair of dirty STL | winding numbers and signed heat method for classification; untangling (Minarcik et al. 2026, as cited in the L13 plan) as an optional repair tool; subgrid marching tets to regenerate | in the L13 plan, plus the new candidate |

### 6.4 Post-processing and images (L21)

- Iso-surfaces and sections of Q-criterion, temperature and VOF fronts on the
  **same tetrahedral decomposition the solver uses** (watertight across
  faces), by marching tetrahedra; the subgrid variant (section 5) for thin
  sheets.
- Streamlines and pathlines integrated with the tracker (exact cell to cell,
  no interpolation to a grid).
- A **headless ray-cast renderer on the D108 hierarchy** for L21 previews
  and for the gallery (no ParaView in the cloud container). Low priority,
  but it comes almost free once ray casting exists; the maintainer's rules on
  images (look at every image, no symbol-like patterns) apply.

## 7. Does each method keep reproducibility and the MPI plan?

| Method | Bitwise (D93) | Multi-rank (MPI plan W3 to W6, SP3) | Float64 |
|---|---|---|---|
| Cell-box BVH | yes (total orders, ties by id) | per-rank tree plus replicated rank-box tree | yes |
| Exact predicates in tracker | yes, and fixes rounding-dependent decisions | both sides of a processor face decide alike | yes |
| Predicate-based ray test, robust traversal | yes | rays stop at the rank boundary, as bundles do | yes |
| Wide BVH, SIMD | result yes (ties); tree shape must be fixed by total orders | unchanged | yes (narrower SIMD) |
| GPU LBVH and refit | yes if the build uses stable sorts and fixed tie-breaks; check against the host build | unchanged | yes |
| Null-collision tracking, majorant grid | yes (draws keyed by bundle and event) | the majorant grid is global (small); a flight crossing a rank boundary needs handover mid-flight: the main risk | yes |
| Stratified emission, Sobol, ratio tracking, MIS | yes (keyed; sums exact or fixed-point) | unchanged, bundles independent of the partition | yes |
| Path guiding | only if learned from fixed-point tallies in fixed order | the learned structure must be rank independent | yes |
| Heat-method fields | linear solves: reproducible to solver tolerance, not bitwise unless the solver is | ordinary scalar solves | yes |
| WoS probes | yes (keyed RNG) | trivially parallel | yes |
| Hardware ray tracing (P9) | **no** (single precision, hardware dependent) | outside the strict mode | no |

## 8. Licences and patents (clean-room check)

Licence files were read from the projects' repositories on 2026-11-24;
nothing else of their code was read.

| Library or algorithm | Licence | Use |
|---|---|---|
| Zombie (Sawhney, WoS, WoSt, geometry queries) | MIT | may be read; write our own |
| geometry-central (Sharp; intrinsic triangulations, heat methods) | MIT | may be read; write our own |
| Embree (Intel, BVH and ray tracing, CPU and SYCL GPU) | Apache-2.0 (licence file and README read); **float32 only** | timing reference only; not a dependency (float64 rule, D60) |
| OpenVDB (sparse level sets, fast sweeping) | Apache-2.0 (the repository notes that the relicensing from MPL-2.0 is complete; check the licence of the exact release before reading) | may be read; write our own |
| pbrt-v4 (the book's renderer) | Apache-2.0 | may be read; the book is the cited source |
| Mitsuba 3 (research renderer) | BSD-style (licence file read) | may be read |
| libigl (fast winding number, Booleans) | MPL-2.0 (licence file read; file-level copyleft) | prefer the paper; do not copy files |
| CGAL (robust geometry, Booleans) | mixed GPL and LGPL parts by package and a commercial licence (the repository file only points to the release licence, not re-read here) | **do not read** (copyleft rule of D201) |
| OptiX (NVIDIA ray tracing) | proprietary SDK licence | not used |
| Star-Engine, EDStar and the Meso-Star libraries (grid-free MC radiation) | not verified here | verify the licence before even reading; the methods are in the published papers |
| p4est and similar adaptive octree libraries | not verified here | verify before reading for L14 |
| Watertight ray-triangle test (Woop, Benthin and Wald 2013) | the paper is open; **patents titled "Watertight ray triangle intersection" exist** (US 10,825,230 and continuations; the claims and the assignee could not be read from this environment, the NVIDIA developer forum thread of the same title suggests NVIDIA) | do not implement this algorithm without a patent check; use the predicate-based test (section 3.3) |
| Ray offsetting for self-intersection (Wächter and Binder 2019) | a published method; a patent application with a similar title exists (not read) | not needed (section 3.3) |

Clean-room consequence: every method above is implemented from the paper
and our own predicates; a licence or patent question stops the item, not
the module.

## 9. Found by search, not read

Titles only, to be read when an item starts; none was used for a ranking:

- TU Delft, a GPU Monte Carlo radiative heat transfer implementation for
  coupling with direct numerical simulation (arXiv:1810.00188).
- "Coupling radiative, conductive and convective heat-transfers in a single
  Monte Carlo algorithm: a general theoretical framework for linear
  situations" (PubMed Central PMC10079137) and "Combined
  conductive-convective-radiative heat transfer in complex geometry using
  the Monte Carlo method" (HAL hal-01882644).
- Wall distance: NASA NTRS 20240006379 (a distance function code paper);
  arXiv:2211.02319 (evaluating a distance function); arXiv:2111.09859
  (high-order differential-equation wall distance). To be read if the
  exact BVH distance ever needs a cross-check beyond the PDE form.
- The Monte Carlo fluid papers of section 5 (University of Waterloo
  preprints).
- Crane's 2025 SIGGRAPH course on grid-free Monte Carlo for PDEs and the CMU
  Monte Carlo Methods course (abstract only): the best single entry point
  when R2 or S4 starts.

## 10. Recommended order

| Step | What | Why this place | Size (agent time) | Needs the maintainer's go |
|---|---|---|---|---|
| 1 | Cell-box BVH, tracker fallback and D52 donor search on it; benchmark against the bin grid on the heat exchanger mesh (D108 part 2a) | sure gain, no dependency, unblocks interpolation, probes and spray injection | one PR, about 1 h code plus CI | yes (code) |
| 2 | Exact predicates (D225), then used in the tracker (3.2) | already planned; makes the tracker exact | per the L13 plan, 2 PRs | yes (code; the L13 go) |
| 3 | Ray casting with the predicate-based test, closest hit and any hit (D108 part 2b) | needed by surface radiation, thin walls, renders; before any radiation | two PRs | yes |
| 4 | Radiation R1 on the walk with stratified emission and Sobol points; **then** the majorant-grid tracker and the MIS estimator as benchmarked options | the largest accuracy-per-cost lever | per the L15 plan; the options not sized yet | yes |
| 5 | L13 benchmark list: add subgrid marching tetrahedra for the penalisation iso-surface and the dirty-STL step | one line in the plan, no cost now | docs only | no |
| 6 | Wide BVH and GPU build (3.4), query ordering | only when ray counts show the need; with P9 | not sized yet | later |
| 7 | Research spikes: path guiding, WoS probes for verification (S4), Hilbert order against RCM | no module is blocked by them | each a small spike | later |
| Not recommended | heat method for the RANS wall distance; Monte Carlo fluids; WoS as P1, flow, film or solid solver; the watertight test of 2013 as published | section 5, 8 | | |

## 11. Questions for the maintainer

1. **First build item.** Start step 1 (cell-box BVH and benchmark) now, or
   wait until the spray or radiation work needs it? Recommendation: start
   with step 1, it is small and measurable and gives the first number for the
   "state of the art" rule.
2. **Null-collision tracker.** Accept that radiation gets *two* estimators
   (walk with exact track length, and null-collision on a majorant grid)
   with the winner chosen by a benchmark on two cases? Recommendation: yes.
3. **Patent rule.** Accept "no algorithm with a known patent family without
   a check, our own predicate-based ray test instead"? Recommendation: yes.
4. **Heat method for the wall distance.** Confirm that D48 stays on the
   exact BVH distance and that the heat method family is judged only inside
   L13 (signed distance)? Recommendation: yes.
5. **Subgrid marching tetrahedra.** Add it to the L13 stage 0 and stage 1
   benchmark lists after reading the full paper? Recommendation: yes.

**Answers (owner 2026-11-24, "go your recommended order"):** all five as
recommended. Question 5 is done in the L13 plan (stage 0 benchmark and
stage 1 iso-surface extraction, #239).

## 12. References

Open literature; each entry was checked against its bibliographic data,
and only the abstract was read for the Crane papers (see the top of the page).

- Attene, M. (2020), Indirect predicates for geometric constructions,
  Computer-Aided Design 126.
- Baktash, H., Gillespie, M., Crane, K. (2025/2026), Subgrid marching
  tetrahedra, ACM Trans. Graph.
- Burley, B. (2020), Practical hash-based Owen scrambling, J. Comput. Graph.
  Techniques 9(4).
- Cherchi, G., Pellacini, F., Attene, M., Scateni, R. (2022), Interactive and
  robust mesh booleans, ACM Trans. Graph. 41(6).
- Chislett, M., Sharp, N., Jacobson, A., Crane, K. (2023), Surface
  simplification using intrinsic error metrics, ACM Trans. Graph. 42(4).
- Corman, E., Crane, K. (2025), Rectangular surface parameterization, ACM
  Trans. Graph.
- Crane, K., Weischedel, C., Wardetzky, M. (2013), Geodesics in heat, ACM
  Trans. Graph. 32(5); (2017), The heat method for distance computation,
  Commun. ACM.
- Dammertz, H., Hanika, J., Keller, A. (2008), Shallow bounding volume
  hierarchies for fast SIMD ray tracing of incoherent rays, Computer
  Graphics Forum 27(4).
- Eisenacher, C., Nichols, G., Selle, A., Burley, B. (2013), Sorted deferred
  shading for production path tracing, Computer Graphics Forum (EGSR) 32(4).
- Eymet, V. et al. (2013), null-collision meshless Monte Carlo, J. Quant.
  Spectrosc. Radiat. Transfer (title shortened).
- Feng, N., Crane, K. (2024), A heat method for generalized signed distance,
  ACM Trans. Graph. 43(4).
- Feng, N., Gillespie, M., Crane, K. (2023), Winding numbers on discrete
  surfaces, ACM Trans. Graph. 42(4).
- Feng, N., Gkioulekas, I., Crane, K. (2026), Points as tori: fast pointwise
  signed distance for point clouds, ACM Trans. Graph. 45.
- Galtier, M. et al. (2013), Integral formulation of null-collision Monte
  Carlo algorithms, JQSRT 125.
- Gillespie, M., Sharp, N., Crane, K. (2021), Integer coordinates for
  intrinsic geometry processing, ACM Trans. Graph. 40(6).
- Gillespie, M., Yang, D., Botsch, M., Crane, K. (2024), Ray tracing harmonic
  functions, ACM Trans. Graph. 43(4).
- Ize, T. (2013), Robust BVH ray traversal, J. Comput. Graph. Techniques 2(2).
- Karras, T. (2012), Maximizing parallelism in the construction of BVHs,
  octrees, and k-d trees, High-Performance Graphics.
- Karras, T., Aila, T. (2013), Fast parallel construction of high-quality
  bounding volume hierarchies, High-Performance Graphics.
- Kutz, P., Habel, R., Li, Y. K., Novak, J. (2017), Spectral and decomposition
  tracking for rendering heterogeneous volumes, ACM Trans. Graph. 36(4).
- Laine, S., Karras, T., Aila, T. (2013), Megakernels considered harmful:
  wavefront path tracing on GPUs, High-Performance Graphics.
- MacDonald, J. D., Booth, K. S. (1990), Heuristics for ray tracing using
  space subdivision, The Visual Computer 6.
- Macpherson, G. B., Nordin, N., Weller, H. G. (2009), Particle tracking in
  unstructured, arbitrary polyhedral meshes for use in CFD and molecular
  dynamics, Commun. Numer. Meth. Engng 25.
- Meister, D., Bittner, J. (2018), Parallel locally-ordered clustering for
  bounding volume hierarchy construction, IEEE Trans. Vis. Comput. Graph. 24(3).
- Meister, D. et al. (2021), A survey on bounding volume hierarchies for ray
  tracing, Computer Graphics Forum 40(2).
- Minarcik, J., Liu, Crane, K., Li (2026), untangling of surface
  self-intersections (as cited in the L13 plan, which lists the authors).
- Bakbouk, G. et al. (2023), mean value caching for walk on spheres (see the
  radiation assessment).
- Miller, B., Sawhney, R., Crane, K., Gkioulekas, I. (2023, 2024, 2025),
  boundary value caching; Walkin' Robin; solving partial differential
  equations in participating media (see the radiation assessment for the full
  entries).
- Möller, T., Trumbore, B. (1997), Fast, minimum storage ray-triangle
  intersection, J. Graphics Tools 2(1).
- Müller, T., Gross, M., Novak, J. (2017), Practical path guiding for
  efficient light-transport simulation, Computer Graphics Forum 36(4).
- Novak, J., Selle, A., Jarosz, W. (2014), Residual ratio tracking for
  estimating attenuation in participating media, ACM Trans. Graph. 33(6).
- Rioux-Lavoie, D. et al. (2022), A Monte Carlo method for fluid simulation,
  ACM Trans. Graph. 41(6).
- Sawhney, R., Crane, K. (2020), Monte Carlo geometry processing, ACM Trans.
  Graph. 39(4); Sawhney et al. (2022, 2023), see the radiation assessment.
- Sharp, N., Crane, K. (2020), A Laplacian for nonmanifold triangle meshes,
  Symp. Geometry Processing.
- Sharp, N., Soliman, Y., Crane, K. (2019), The vector heat method, ACM Trans.
  Graph. 38(3); Navigating intrinsic triangulations, ACM Trans. Graph. 38(4).
- Shewchuk, J. R. (1997), Adaptive precision floating-point arithmetic and
  fast robust geometric predicates, Discrete and Computational Geometry 18.
- Szirmay-Kalos, L., Toth, B., Magdics, M. (2011), Free path sampling in high
  resolution inhomogeneous participating media, Computer Graphics Forum 30(1).
- Veach, E., Guibas, L. J. (1995), Optimally combining sampling techniques for
  Monte Carlo rendering, SIGGRAPH.
- Vose, M. D. (1991), A linear algorithm for generating random numbers with a
  given distribution, IEEE Trans. Software Eng. 17(9); Walker, A. J. (1977),
  An efficient method for generating discrete random variables with general
  distributions, ACM Trans. Math. Software 3(3).
- Wächter, C., Binder, N. (2019), A fast and robust method for avoiding
  self-intersection, in Ray Tracing Gems.
- Wald, I. (2007), On fast construction of SAH-based bounding volume
  hierarchies, IEEE Symp. Interactive Ray Tracing.
- Wilkie, A., Nawaz, S., Droske, M., Weidlich, A., Hanika, J. (2014), Hero
  wavelength spectral sampling, Computer Graphics Forum 33(4).
- Woodcock, E. et al. (1965), Techniques used in the GEM code for Monte Carlo
  neutronics calculations, ANL-7050.
- Woop, S., Benthin, C., Wald, I. (2013), Watertight ray/triangle
  intersection, J. Comput. Graph. Techniques 2(1) (see the patent note).
- Ylitie, H., Karras, T., Laine, S. (2017), Efficient incoherent ray traversal
  on GPUs through compressed wide BVHs, High-Performance Graphics.
- Yoon, S.-E., Lindstrom, P., Pascucci, V., Manocha, D. (2005), Cache-oblivious
  mesh layouts, ACM Trans. Graph. 24(3).
