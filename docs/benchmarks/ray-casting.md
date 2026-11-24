# Ray casting: exact segment–triangle test against Möller–Trumbore

Benchmark of 2026-11-24 for D108 part 2 (pull request #RAY), step 3 of the
[graphics-algorithms assessment](../assessments/graphics-algorithms-for-modules.md#33-ray-casting-on-the-bvh-and-what-to-build-it-from).

**What is compared.** The test that decides whether a ray or a segment
meets a triangle of a wall surface:

- **exact** (ours, `foehncraft_mesh::query::segment_triangle`): the
  segment from `p` to `q` crosses the plane when the two signed volumes of
  `p` and `q` with the triangle differ in sign, and its line passes the
  closed triangle when the three volumes of the line with the triangle
  edges have no two of opposite sign (O'Rourke 1998, §7.3; Segura & Feito
  2001). All five volumes are the exact orientation predicate of D225
  (Shewchuk 1997): a floating-point filter, exact arithmetic only near
  zero. Two triangles with a common edge evaluate its volume as exact
  negations, so no ray slips between them;
- **Möller–Trumbore** (1997), the common floating-point test of renderers,
  written into the benchmark from the paper: barycentric weights from cross
  products, accepted when both are in [0, 1] (no tolerance) or in
  [−10⁻⁹, 1 + 10⁻⁹] (tolerance, which accepts a ray in both neighbours of
  an edge).

The watertight test of Woop, Benthin and Wald (2013) is not compared: it
carries patents (assessment section 8).

**How to run it.**

```bash
cargo run --release -p foehncraft-mesh --example ray_bench
```

The surface is the closed boundary of a generated box mesh (node
perturbation 0.2, polygon faces split into fans around the node average),
turned about an oblique axis so that no wall lies in a coordinate plane.

- **Leaks**: rays from four interior cell centroids aimed exactly at every
  surface corner (mesh nodes and fan centres) and at two random points on
  every surface edge. The box is convex, so each ray leaves it exactly at
  the aimed corner or edge and only the triangles there can catch it; a
  leak is a ray that none of them catches.
- **Test cost**: nanoseconds per test on 2 million random ray–triangle
  pairs (mostly misses, as inside a traversal).
- **First hit**: microseconds per `Bvh::ray_first_hit` from a random cell
  centroid in a random direction (100 000 rays, all must hit), and per
  `Bvh::first_hit_brute_force` (200 rays, the same hits as the hierarchy).

Single thread, 4-core cloud container (Linux x86-64), release build.

| Mesh | Triangles | Rays | Leaks exact | Leaks Möller–Trumbore | Leaks Möller–Trumbore 1e-9 | Exact test (ns) | Möller–Trumbore (ns) | First hit on the hierarchy (µs/ray) | Brute force (µs/ray) |
|---|---|---|---|---|---|---|---|---|---|
| hex 64³ | 98304 | 1376264 | 0 | 182584 | 0 | 77.3 | 37.9 | 3.55 | 1914 |
| tet 24³ | 6912 | 96776 | 0 | 12911 | 0 | 26.2 | 19.9 | 2.27 | 202 |
| poly 20³ | 29760 | 416648 | 0 | 39152 | 0 | 58.9 | 32.4 | 2.41 | 479 |

**What it shows.**

- The exact test never leaks. Möller–Trumbore without a tolerance loses
  9 to 13 % of the rays aimed at corners and edges: rounding of the two
  barycentric weights is independent in the two neighbours, so both can
  reject. Random rays hit an edge or a corner rarely, so a renderer sees
  this as an occasional speckle; a Monte Carlo radiation tally or a parcel
  tracker sees it as energy or mass that leaves through a closed wall.
- A tolerance stops the leaks but makes a ray hit both neighbours of an
  edge (double counting for parity and for crossing counts) and depends on
  the triangle size. The exact test has no tolerance to choose.
- The exact test costs 1.3 to 2 times a Möller–Trumbore test. A first hit
  on the hierarchy takes 2 to 4 µs per ray on 7 000 to 100 000 triangles,
  540 times faster than testing every triangle on the largest surface.
- Rays started on a wall need no offset: the segment is half open at its
  start, so a start point in the plane of a triangle does not meet it
  (tested in `crates/foehncraft-mesh/tests/ray.rs`).

**Possible later speed-up.** The three edge volumes share the line through
`p` and `q`. In Plücker coordinates of that line (Amanatides & Choi 1995;
Jones 2000) each volume is one inner product of six numbers; with a static
error bound as the filter and the exact predicate as the fallback the test
would keep its decisions and cost about as much as Möller–Trumbore. Not
done until a user (spray wall impact, surface-to-surface radiation)
measures the test as the bottleneck.
