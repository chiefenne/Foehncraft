# Point search: hierarchy over cell boxes against the uniform grid

Benchmark of 2026-11-24 for D108 part 2 (pull request #235), step 1 of the
[graphics-algorithms assessment](../assessments/graphics-algorithms-for-modules.md#31-hierarchy-over-cell-boxes-replaces-the-uniform-bin-fallback).

**What is compared.** The search for the cell that holds a point when no
hint cell is known (probes, line samples, the pressure reference cell,
injection points, the fallback of the particle walk of D47):

- **grid** (before): a uniform grid of cell bounding boxes, bin size twice
  the mean cell size, at most 1024 bins per axis; the cells of the bin of
  the point are tested in index order (`Locator::locate_search_grid`);
- **tree** (now): a bounding volume hierarchy over the cell boxes (binned
  surface-area heuristic, at most four cells per leaf, the build of the
  wall-distance hierarchy), every cell whose box holds the point is tested
  and the smallest index that contains it is returned
  (`Locator::locate_search`).

**How to run it.**

```bash
cargo run --release -p foehncraft-mesh --example locate_bench -- \
  examples/meshes/naca0012.su2 examples/meshes/diffuser_coarse.su2 \
  examples/meshes/ubend_coarse.su2 examples/meshes/tjunction_gci_coarse.msh
```

Queries: every cell centroid (at most 200 000) and as many random points in
the bounding box of the mesh (outside points included for non-convex
domains). Generated meshes have a node perturbation of 0.1; "graded" is a
symmetric grading of the y spacing with the given ratio of the middle to
the end intervals.

## Results (cloud container, Intel Xeon at 2.1 GHz, 4 cores)

Queries run on one thread; the tree build uses the 4 cores (the grid build
is serial).

| Mesh | Cells | Grid build (s) | Tree build (s) | Queries | Grid (µs per query) | Tree (µs per query) | Speed-up | Candidate cells per query (tree) | Different answers |
|---|---|---|---|---|---|---|---|---|---|
| hex 64³ uniform | 262 144 | 0.116 | 0.097 | 400 000 | 12.95 | 1.91 | 6.8 | 1.2 | 0 |
| hex 64³ graded 1e4 in y | 262 144 | 0.095 | 0.102 | 400 000 | 26.72 | 1.33 | 20.0 | 1.2 | 0 (148 before #FC, see below) |
| tet 24³ graded 1e3 in y | 82 944 | 0.032 | 0.035 | 165 888 | 8.35 | 1.07 | 7.8 | 6.5 | 0 |
| poly 20³ graded 1e3 in y | 9 261 | 0.011 | 0.026 | 18 522 | 50.78 | 4.98 | 10.2 | 2.3 | 0 |
| `naca0012.su2` (2D airfoil, boundary layer) | 6 420 | 0.002 | 0.005 | 12 840 | 6.95 | 0.87 | 8.0 | 1.6 | 0 |
| `diffuser_coarse.su2` | 6 960 | 0.002 | 0.004 | 13 920 | 2.66 | 0.38 | 7.0 | 1.1 | 0 |
| `ubend_coarse.su2` | 8 832 | 0.002 | 0.006 | 17 664 | 3.08 | 0.35 | 8.7 | 1.3 | 0 |
| `tjunction_gci_coarse.msh` (tetrahedra) | 20 483 | 0.007 | 0.017 | 40 966 | 9.31 | 1.72 | 5.4 | 4.4 | 0 |

Build on a larger mesh (hexahedra 128³, 2 097 152 cells): grid 0.65 s, tree
0.91 s (the serial tree build took 2.9 s; the parallel build splits large
ranges into subtrees built at the same time, and bins with exact chunked
counts and min/max boxes).

**Reading.** The tree is 5 to 20 times faster per query; the gain grows
with the grading, because the grid sizes its bins by the mean cell and
then tests many small cells per bin. The tree tests 1 to 7 candidate cells
per query on every mesh. Its build costs about as much as the grid's on the
large meshes and two to three times as much on the small ones (milliseconds),
paid once per mesh on the first query.

**The 148 different answers** (0.04 % of the queries, all on the generated
hexahedra with a grading of 1e4 and node perturbation) are cells that the
containment test accepts although the point lies outside the cell's box: on
these ribbon-thin, twisted side faces (thickness 1e-5 against a perturbation
of 1e-3 along the face) the area-weighted face centroid of
`geometry::polygon_area_centroid` falls far outside the face (z = 0.19 for a
face between z = 0 and 0.017), because the fan triangles' projected areas
change sign and nearly cancel. The cell apex, and with it the tetrahedral
decomposition, then leaves the cell. The grid tests such a cell when it
shares the point's bin; the tree never tests a cell whose box does not hold
the point, so its answer is the right one. It does not occur on the
meshes of the repository (0 differences).

**Fixed in #FC (2026-11-24).** When the fan triangles of a face turn
against each other (a projected area of zero or less) and the area-weighted
centroid leaves the bounding box of the face nodes, the centroid now uses
the triangle area magnitudes as weights, which keeps it on the face. Planar
faces, convex or not, keep the exact area centroid, so no other mesh
changes (the golden geometry values are unchanged). With the fix the grid,
the tree and the brute-force search agree on every query of every mesh
above (rerun of `locate_bench`: 0 different answers, also on the graded
hexahedra).

## Determinism and limits

- The tree build orders with total orders only and the answer is the
  smallest cell index that contains the point, so the result is the same on
  every machine and thread count (D93); a unit test checks that two builds
  give the same tree and that the tree, the grid and a brute-force search
  agree on uniform and graded meshes of every cell type.
- The build of the wall-distance hierarchy shares this code: it splits
  with a stable partition instead of a sort (linear instead of n log n per
  level, the same split as before) and builds large subtrees in parallel;
  its queries give the same answers (the wall-distance tests are
  unchanged).
- The build is parallel above 16 384 items; where it splits into parallel
  subtrees depends only on the range sizes, so the tree is the same for any
  number of threads (a unit test builds it on one thread and on all and
  compares node by node). The queries run on the caller's thread; a
  parallel query loop is trivial (the tree is read-only) when a user needs
  it.
