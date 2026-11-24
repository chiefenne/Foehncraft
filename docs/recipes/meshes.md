# Meshes

Part of the [Recipes](../recipes.md): generate, read (Gmsh, SU2, CGNS, OpenFOAM, Fluent .msh, VTK), convert and partition meshes, update the mesher's golden values, run flows on poor Gmsh meshes.

| # | Recipe |
|---|---|
| 7 | [Generate a test mesh and check its quality](#mesh-generate) |
| 8 | [Read a Gmsh mesh and convert it into polyhedra](#mesh-gmsh) |
| 9 | [Partition a mesh and test the halo exchange](#mesh-partition) |
| 11 | [Update the golden values of the mesher](#golden) |
| 36 | [Run a flow on raw Gmsh tetrahedra (skewed faces)](#skewed-tets) |
| 127 | [Read an SU2 mesh and map its markers to boundary conditions](#mesh-su2) |
| 129 | [Read a CGNS or OpenFOAM mesh](#mesh-cgns-openfoam) |
| 134 | [Read a Fluent `.msh` or VTK mesh, write a VTU](#mesh-msh-vtk) |
| 138 | [Run a rotor/stator sliding case from a multi-zone SU2 or CGNS file](#mesh-multizone) |
| 189 | [Check which cells the velocity gradient and Hessian are exact on](#operator-check) |
| 200 | [Check that every rank has the bits of one rank (canonical local mesh)](#canonical-ranks) |
| 219 | [Sum over ranks exactly (and check it)](#exact-reductions) |
| 230 | [Solve a scalar on ranks (and check it against one rank)](#scalar-ranks) |
| 231 | [Use the multigrid preconditioner on ranks](#multigrid-ranks) |
| 234 | [Solve a flow on ranks (and check it against one rank)](#flow-ranks) |
| 244 | [Check a turbulent, heated or transient flow on ranks](#flow-models-ranks) |

---

<a id="mesh-generate"></a>
## 7. Generate a test mesh and check its quality

The built-in mesher (D34) makes a box of hexahedra, tetrahedra, polyhedra
or mixed cells, in 3D or 2D (one cell thick), graded, refined and
perturbed. Keys `[mesh.generate]` ([Case catalogue](../cases/meshes.md#mesh-case-file-keys)).

```
cargo build --release -p foehncraft-cli
target/release/foehn mesh examples/cases/mesh/poly_box.toml
cat output/poly_box/mesh_report.md
```

The report gives the verdict (valid or not), the cell types, the closure
|Σ S| / Σ |S| of the worst cell, the total volume against the volume
enclosed by the boundary, the boundary areas, non-orthogonality, skewness,
face warp, aspect and volume ratio with the cell or face where the maximum
sits, and the point-location test. `mesh.vtkhdf` holds the mesh with these
fields per cell; open it in ParaView (5.13 or newer reads VTKHDF 2.4 with
polyhedra) or render a picture:

```
pip install vtk
python3 -m tools.mesh.render_mesh output/poly_box/mesh.vtkhdf poly.png --cut --field non_orthogonality
```

A grid-convergence family: the same case with `refine = 0, 1, 2` gives
nested meshes (every interval halved), the grading kept.

[↑ back to top](#meshes)

---

<a id="mesh-gmsh"></a>
## 8. Read a Gmsh mesh and convert it into polyhedra

Needs Gmsh 4.x (`pip install gmsh` or gmsh.info). The reader takes MSH 4.1,
ASCII or binary, first-order elements (tetrahedra, hexahedra, prisms,
pyramids), physical groups and periodic links (D33).

```
gmsh -3 examples/meshes/tjunction.geo -format msh41 -o examples/meshes/tjunction.msh
target/release/foehn mesh examples/meshes/tjunction.msh             # defaults, no case file
target/release/foehn mesh examples/cases/mesh/tjunction_tet.toml    # RCM, 4 parts
target/release/foehn mesh examples/cases/mesh/tjunction_poly.toml   # converted into polyhedra
```

In your own case file: `mesh.file` (relative to the case file),
`mesh.scale = 0.001` for a mesh in millimetres, `mesh.polyhedral = true`
for the conversion of a tetrahedral mesh into polyhedra (one cell per node,
same volume), `mesh.renumber = "rcm"` for a band-limited cell order, and
`[mesh.boundaries]` to rename physical groups. Give every boundary surface a
physical group in Gmsh; surfaces without one become `surface_<tag>`.

[↑ back to top](#meshes)

---

<a id="mesh-partition"></a>
## 9. Partition a mesh and test the halo exchange

```
target/release/foehn mesh examples/cases/mesh/poly_partition.toml
target/release/foehn mesh examples/cases/mesh/poly_partition.toml --parts 8 --method rcb
```

Keys `partition.parts`, `partition.method` (`metis` or `rcb`),
`partition.halo_layers`. The report lists cells, halo cells, processor
faces and neighbours per rank, the imbalance and the cut faces;
`mesh.vtkhdf` has the field `rank`.

The halo exchange over MPI (needs MPI as in recipe 4): every rank builds
its part with a two-layer halo and fills it from its neighbours.

```
cargo build --release -p foehncraft-mesh --features mpi --example mpi_halo
mpirun --oversubscribe -n 4 target/release/examples/mpi_halo
```

Prints `MPI halo exchange: PASS` when the halo values, volumes and centroids
equal the global mesh and the volume sum gives the same bits on every rank.

[↑ back to top](#meshes)

---

<a id="golden"></a>
## 11. Update the golden values of the mesher

`crates/foehncraft-mesh/tests/golden/mesher.json` holds the cell, face and
node counts, volumes, quality metrics and partition statistics of every
combination of mesh type, dimension, perturbation, refinement and 1, 2 and
4 ranks (R4). A difference fails the test. Only for an intended change:

```
FOEHNCRAFT_BLESS=1 cargo test -p foehncraft-mesh --test golden
git diff crates/foehncraft-mesh/tests/golden/mesher.json
```

and say in the commit message why the values changed.

[↑ back to top](#meshes)

---

<a id="skewed-tets"></a>
## 36. Run a flow on raw Gmsh tetrahedra (skewed faces)

Gmsh's tetrahedra without its optimisation pass have faces of a skewness
above 0.9, mostly at curved walls; the steady solver diverges on them with
the least-squares pressure gradient ([gallery 26](../gallery/meshes.md#skewed-tets)).
First choice: optimise the mesh. Add `Mesh.OptimizeNetgen = 1;` to the .geo
file (or `-optimize_netgen` on the command line) and check the count in the
quality report:

```
gmsh -3 examples/meshes/tjunction_netgen.geo -format msh41 -o examples/meshes/tjunction_netgen.msh
target/release/foehn mesh examples/meshes/tjunction_netgen.msh
grep "skewness above 0.9" output/tjunction_netgen/mesh_report.md
```

When the mesh has such faces, the flow solver treats them itself
(`solver.remediation = "auto"`, the default): it switches to the Gauss
pressure gradient, the Venkatakrishnan limiter and the bounded face
velocity of the mass flux (`solver.face_velocity = "bounded"`) for the
whole run and
says so at the start of `flow.log` and in `flow.md`, `flow.json`
(`remediation` of each solve) and the diagnostics (`MESH001`, `MESH002`):

```
target/release/foehn flow examples/cases/flow/pipe_gmsh_raw.toml
grep -E "mesh quality|remediation" output/pipe_gmsh_raw/flow.log
```

Keys (`[solver]`): `remediation = "auto"` or `"off"` (keep the schemes as
set; the run on `pipe_gmsh_raw` then diverges and `divergence.md` names the
skewed faces), `skewness_limit = 0.9` (a smaller value treats milder
faces), `pressure_gradient = "least-squares"` or `"gauss"` (set by hand,
also on a good mesh), `face_velocity = "interpolated"` or `"bounded"`.
A limiter set in the case stays, and first-order `upwind` needs none (the
bounded face velocity then takes the Venkatakrishnan limiter for itself).
A turbulent run that still grows k in a few cells of a poor mesh can also
take `transposed_stress = "trace-free"` (not switched automatically).
The test case of the bounded face velocity is a coarse stirred tank on raw
tetrahedra (`pip install gmsh`, then `python3 examples/studies/skewed_tets.py`). Least squares is more accurate on good meshes (the
cylinder lift at Re = 20 is 3.5 % lower with the Gauss gradient), so treat
the remediation as a remedy and fix the mesh for production runs.

```
target/release/foehn flow examples/cases/flow/pipe_gmsh_raw.toml --set solver.remediation=off \
    --out output/pipe_gmsh_raw_off
cat output/pipe_gmsh_raw_off/divergence.md
```

[↑ back to top](#meshes)

---

<a id="mesh-su2"></a>
## 127. Read an SU2 mesh and map its markers to boundary conditions

The SU2 reader (D187) takes native `.su2` files, ASCII, one zone, 2D or
3D, first-order elements: tetrahedra, hexahedra, prisms and pyramids in 3D
with triangle and quadrilateral markers; triangles and quadrilaterals in 2D
with line markers. A 2D file is extruded one cell layer in z (triangles to
prisms, quadrilaterals to hexahedra), and the planes become the boundaries
`zmin` and `zmax`, as for the built-in 2D meshes. Every `MARKER_TAG` is a
boundary of that name.

```
target/release/foehn mesh examples/meshes/naca0012.su2               # defaults, no case file
target/release/foehn mesh examples/cases/mesh/naca0012_su2.toml      # thickness, RCM, 4 parts
target/release/foehn flow examples/cases/flow/naca0012_su2_re1000.toml
```

(`examples/meshes/naca0012.su2` is a PyAero mesh in the repository, see
[gallery 127](../gallery/meshes.md#su2-airfoil).) In your own case file:

```toml
[mesh]
file = "wing.su2"          # the format follows from the extension, or:
format = "su2"             # gmsh | su2 | fluent_msh | cgns | vtk | openfoam, for a file with another extension
thickness = 0.01           # 2D files: layer thickness, in the file's units
scale = 0.001              # millimetres to metres (after the extrusion)

[mesh.boundaries]          # marker names to the names of the conditions
farfield = "inlet"
wing = "wall"

[boundary.wall]
type = "wall"
```

Without `thickness` the layer is as thick as the median edge of the 2D
cells. Boundary faces that no marker names go to `default` with a warning.
A marker called `zmin` or `zmax` in a 2D file merges with the extruded
plane of that name: rename it with `[mesh.boundaries]`. Multi-zone files
(`NZONE` > 1) are refused: write one file per zone. The format layer
picks the reader from `mesh.format`, else the extension (`.msh`, `.su2`),
else the first lines of the file; a new format is one reader in
`crates/foehncraft-mesh/src/format/` (module docs of `format`).

[↑ back to top](#meshes)

---

<a id="mesh-cgns-openfoam"></a>
## 129. Read a CGNS or OpenFOAM mesh

Two more readers sit behind the same interface (D188).

- **OpenFOAM.** Give the case directory, its `constant` directory or the
  `polyMesh` directory itself: `file = "channel"`. The files `points`,
  `faces`, `owner`, `neighbour` and, if present, `cellZones` are read, ASCII
  or binary (32 or 64 bit labels, either byte order), with polyhedral cells
  kept as they are. Patches become boundaries by name, `cellZones` become
  cell zones. `empty` patches stay ordinary boundaries (the reader says so;
  use a `symmetry` boundary for 2D work). Decomposed cases (`processor`
  patches) and `.gz` files are refused: reconstruct and `gunzip` first.
- **CGNS.** HDF5 files (CGNS 3 or 4 as written by Gmsh and others), one
  base, one unstructured zone, 3D cells TETRA_4, PYRA_5, PENTA_6, HEXA_8
  and polyhedra as NGON_n/NFACE_n. The BCs of `ZoneBC`, given on faces,
  become boundaries; for Gmsh files the physical-group name is taken from
  the family. Old ADF files are refused: convert them with `adf2hdf`.

```
target/release/foehn mesh crates/foehncraft-mesh/tests/data/openfoam/dual_polyhedra
target/release/foehn mesh crates/foehncraft-mesh/tests/data/cgns/box_tet.cgns
```

```toml
[mesh]
file = "channel"            # an OpenFOAM case, or "wing.cgns"

[mesh.boundaries]           # patch or BC names to the names of the conditions
walls = "wall"
```

Faces that no patch or BC names go to `default` with a warning. Multi-zone
and 2D CGNS files and structured zones are refused.

[↑ back to top](#meshes)

---

<a id="mesh-msh-vtk"></a>
## 134. Read a Fluent `.msh` or VTK mesh, write a VTU

Two more readers (D191, D192) and a writer.

- **Fluent `.msh`.** The ASCII, 3D, face-based text file that many
  commercial meshers and solvers write: nodes, cells and faces with the two
  cells on each side, zones with names. Cell shapes do not matter, so
  polyhedral meshes read as they are. Every face zone that is not of the type
  `interior` is a boundary of its name, the cell zones are the cell zones.
  The extension `.msh` is shared with Gmsh: the first bytes tell them apart,
  and `format = "fluent_msh"` forces the reader. Binary files and 2D files
  are refused with the hint what to do.
- **VTK.** Legacy `.vtk` (ASCII, binary) and XML `.vtu` (ASCII, base64,
  appended, zlib-compressed: what ParaView saves). The volume cells are
  tetrahedra, hexahedra, wedges, pyramids and polyhedra; triangles and
  quadrilaterals are boundary faces. A grid has numbers where we have
  names: the cell-data array `region` (or `gmsh:physical`, `CellEntityIds`,
  `PatchID`, ...) numbers the zones (volume cells) and the boundaries
  (surface cells); name them with the field-data string arrays
  `zone_names` and `boundary_names` (entries `3=inlet`), else a boundary is
  `boundary3`. A file with surface cells that carry no region number puts
  them in `default`.
- **Large VTK files.** The reader streams the file: it never holds the file
  or the text of an array whole, and it keeps node ids as 32-bit integers
  and the faces of the polyhedra in one flat array. A mesh of 1 M
  polyhedra (14 faces each, as KoflCraft writes them) peaks at about 1 GB
  while it is read, whatever the encoding, so 10 M cells need about 10 GB.
  Appended raw binary is the fastest to read (1 M cells in 5 s, ASCII in
  8 s); prefer it for meshes of millions of cells.
- **Writing a VTU.** `output.vtu = true` in the mesh case, or `--vtu`:

```
target/release/foehn mesh crates/foehncraft-mesh/tests/data/fluent_msh/dual_polyhedra.msh --vtu
target/release/foehn mesh crates/foehncraft-mesh/tests/data/vtk/channel_mixed_zlib.vtu
```

The `mesh.vtu` of the first command holds the 81 polyhedra and the boundaries
`xmin` to `zmax` with their names; ParaView, and `foehn` itself, read it back.

```toml
[mesh]
file = "channel.msh"        # Fluent .msh or Gmsh, by the content
# format = "fluent_msh"   # to force it

[mesh.boundaries]           # zone or region names to the names of the conditions
boundary3 = "inlet"         # a VTK region number without a name
```

[↑ back to top](#meshes)

<a id="mesh-multizone"></a>
## 138. Run a rotor/stator sliding case from a multi-zone SU2 or CGNS file

A machine mesh has one zone for the rotor and one for the stator, and the
two do not share nodes. SU2 (`NZONE= 2`, one `IZONE=` block per zone) and CGNS
(several `Zone_t`) store it that way, and both are read now (D195):

- all zones come into **one mesh**, every zone a cell zone: `zone1`, `zone2`
  for SU2, the `Zone_t` names for CGNS (a single zone is still `fluid`);
- the nodes of different zones **stay apart**, so the interface is two
  ordinary boundaries, one per side: the markers (SU2) or BCs (CGNS) that
  the mesher put there, e.g. `rotor_interface` in the rotor zone and
  `stator_interface` in the stator zone (rename with `[mesh.boundaries]` if
  the mesher called them otherwise). That is the mesh `[[sliding]]` wants:
  two boundaries, no interior face between the zones, and the two sides may
  have different pitch and node positions;
- markers or BCs of the same name in several zones (`wall`) are one boundary;
- the interface connectivity entries that some writers add are not read.

```toml
[mesh]
file = "rotor_stator.su2"   # NZONE= 2: zone1 the rotor, zone2 the stator

[[sliding]]
zone = "zone1"              # the cell zone that turns
axis = [0.0, 0.0, 1.0]
omega = 10.0
rotor_interface = "rotor_interface"    # the two interface boundaries of the file
stator_interface = "stator_interface"
```

`mesh.merge_nodes = true` joins nodes that coincide (within 1e-8 of the
bounding-box diagonal). Use it where zones **touch conformally and do not move
against each other**: the interface faces become interior faces and the zones
are one connected mesh (the two interface boundaries then lie inside and are
ignored, with a warning). Do not use it for the zones of an imprinted sliding run (the
two interface boundaries have to stay); it is for conformal stepping on a
conformal interface.

Multi-zone SU2 files are 3D only. Try the case of gallery 138:

```
python3 -m tools.mesh.taylor_couette_zones examples/meshes/taylor_couette_zones.su2
target/release/foehn flow examples/cases/flow/taylor_couette_su2_zones.toml
```

[↑ back to top](#meshes)

---

<a id="operator-check"></a>
## 189. Check which cells the velocity gradient and Hessian are exact on

The curvature correction of SST and the QCR differentiate the velocity
gradient once more; they converge with the mesh only where the gradient is
exact for a quadratic field (D236). Every run checks this before the solve;
an agent or a script can ask before running:

```
target/release/foehn flow examples/cases/flow/cavity_sst_operator_check.toml --check --json
```

The report has `operators` with `linear`, `gradient` and `hessian`, each
with `max`, `mean`, `flagged` (the cells above `tolerance`) and `worst`
(up to 50 cells with index, error and centroid), and the diagnostics
`MESH005` (with `turbulence.curvature_correction`) and `MESH006` (with
`turbulence.qcr`); `MESH004` on any case means a broken stencil. In a run
the same lines go to the log and `flow.md`, `LIM008` to `result.json`, and
`operators.csv` holds `cell,x,y,z,linear,gradient,hessian` for every cell
(open it in ParaView with the cell centroids, or join it with
`fields.vtkhdf` by the cell index). The tolerance is
`solver.operator_tolerance` (0.01):

```
target/release/foehn flow examples/cases/flow/cavity_sst_operator_check.toml \
    --set solver.operator_tolerance=0.05
```

What to do with flagged cells: regular hexahedra or polyhedra in the core
of the flow, or (the next steps of D236) the quadratic least-squares
operators on the flagged cells and a local remesh of the cells still
flagged; converting tetrahedra to polyhedra (D54) is the other way.

[↑ back to top](#meshes)

---

<a id="canonical-ranks"></a>
## 200. Check that every rank has the bits of one rank (canonical local mesh)

Every partitioned `foehn mesh` run checks the local meshes against the
one-rank mesh bit for bit (MPI step W1, D93 (1), (3)) and stops with an
error if any face or cell differs:

```
target/release/foehn mesh examples/cases/mesh/canonical_ranks.toml
target/release/foehn mesh examples/cases/mesh/canonical_ranks.toml --parts 4 --method rcb
```

`mesh_report.md` has the columns `flipped faces` (0 since MPI step W3:
faces keep their global orientation) and `faces, cells differing from one
rank` (always `0, 0`) in the partition table;
`mesh_report.json` has the list `canonical` with `rank`, `flipped_faces`,
`differing_faces` and `differing_cells` for an agent or a script;
`mesh.vtkhdf` has the cell fields `rank` and `flipped_faces`.

The same check across MPI processes, with the rows of a diffusion matrix
(needs MPI as in recipe 4; `SPIKE_KIND` is `poly`, `hex`, `tet` or `mixed`):

```
cargo build --release -p foehncraft-mesh --features mpi --example mpi_matrix
for n in 1 2 3 4; do SPIKE_KIND=tet mpirun --oversubscribe -n $n target/release/examples/mpi_matrix; done
```

Prints `BITWISE EQUAL TO ONE RANK` and exits with 1 otherwise. In code:
`decompose::canonical_check(&global, &Geometry::compute(&global), &locals)`;
the reference is the mesh handed to `decompose` (after any renumbering), so
runs without a decomposition do not change.

[↑ back to top](#meshes)

---

<a id="exact-reductions"></a>
## 219. Sum over ranks exactly (and check it)

Every partitioned `foehn mesh` run also sums the domain volume, its first
moment and the boundary area per rank into exact accumulators, merges them
as an MPI all-reduce would, and stops with an error if any result differs
in a bit from the one-rank sum (MPI step W2, D31 (1), D93 (4)):

```
target/release/foehn mesh examples/cases/mesh/exact_reductions.toml
target/release/foehn mesh examples/cases/mesh/exact_reductions.toml --parts 13 --method rcb
```

`mesh_report.md` has the table "Exact reductions over the ranks" with the
exact value, `bitwise = one rank` (always yes) and the old rank-order
pairwise sum with its distance in ulps; `mesh_report.json` has the list
`reductions` (`name`, `ranks`, `one_rank`, `bitwise_equal`,
`pairwise_ranks`, `pairwise_ulps`).

Across MPI processes, with a dot product, a residual norm with the halo from
the exchange, a boundary flux, min and max and a force in one all-reduce
(needs MPI as in recipe 4):

```
cargo build --release -p foehncraft-mesh --features mpi --example mpi_reductions
for n in 1 2 3 4; do mpirun --oversubscribe -n $n target/release/examples/mpi_reductions; done
```

Prints `BITWISE EQUAL TO ONE RANK` and exits with 1 otherwise. In code:

```rust
use foehncraft_core::exact::{ExactSum, ExactSumIter as _};
let mut acc = ExactSum::new();
acc.add_products(&r, &r);                    // or add_slice, add
let norm = comm.sum_exact(&acc).sqrt();      // the same bits on any rank count
let area: f64 = faces.iter().map(|f| f.area).exact_sum();   // one rank
let [fx, fy] = comm.sum_exact_many(&[ax, ay])[..] else { unreachable!() };
```

The cost per value of the exact dot product against the old block sum:
`cargo run --release -p foehncraft-parallel --example exact_reduction_cost`.

[↑ back to top](#meshes)

---

<a id="scalar-ranks"></a>
## 230. Solve a scalar on ranks (and check it against one rank)

List rank counts under `[run]` of a scalar case: every variant is solved
once on all cells and then on each rank count, with the ranks as threads of
this process (recursive coordinate bisection, one halo layer), and the run
fails unless φ is bitwise the first backend's (MPI step W3, D93):

```toml
[run]
ranks = [2, 3, 5]
```

```
target/release/foehn scalar examples/cases/scalar/ranks_oblique_step_2d.toml
```

`scalar.md` and `scalar.json` list a solve per rank count (`backend` =
`cpu, N ranks`, `bitwise_equal`); `fields.vtkhdf` gets the cell field
`rank` of the largest count. Since MPI step W4 the multigrid preconditioner
runs on ranks too (recipe 231).

Across MPI processes, the solution, gradient, limiter, every matrix entry
and the residual history of five cases (polyhedra, mixed cells,
tetrahedra, hexahedra; linear upwind, bounded central, BiCGStab, CG, CG +
multigrid, the reduction multigrid) compared bit for bit with one rank (needs MPI as in recipe 4):

```
cargo build --release -p foehncraft-fv --features mpi --example mpi_scalar
for n in 1 2 3 4; do mpirun --oversubscribe -n $n target/release/examples/mpi_scalar; done
```

In code, one rank's solve on its part and the check
(`foehncraft_fv::ranks`; `ThreadComm::run` gives ranks as threads in tests):

```rust
let local = extract(&global, &part, comm.rank(), 1)?;      // faces in global orientation
let geo = local.geometry(&*comm);                           // halo volumes and centroids
let d = Discretisation::new(&local.mesh, &geo, &problem)?;
let mut be = CpuBackend::on_ranks(d, threads, Ranks::new(comm, &local, &part))?;
be.upload(V::Phi, &phi0)?;                                  // own cells
let outcome = solver::solve(&mut be, &settings, |_| {})?;  // the one-rank solver, unchanged
```

[↑ back to top](#meshes)

---

<a id="multigrid-ranks"></a>
## 231. Use the multigrid preconditioner on ranks

Nothing changes in the case file: `preconditioner = "amg"` (CG on a
symmetric matrix) runs on ranks as on one rank, and `[run] ranks` checks it
bit for bit (MPI step W4, D93):

```toml
[solver]
linear_solver = "cg"
preconditioner = "amg"
amg = { smoother = "chebyshev", degree = 2, precision = "mixed" }

[run]
ranks = [2, 3, 5]
```

```
target/release/foehn scalar examples/cases/scalar/ranks_poisson_amg_poly_box.toml
```

The setup is replicated: every rank gathers the matrix in global order and
builds the whole hierarchy, so the aggregates are those of one rank and the
setup time does not fall with the rank count (`multigrid setup` in the
phase times of `scalar.md`). The cycle is distributed: each rank smooths,
restricts and prolongs its own rows, exchanges the ghosts of each level and
gathers the coarsest level (at most `max_coarse` unknowns) to solve it
itself. The reduction multigrid of the coupled momentum sweep (D179) uses
the same plans.

In code (`foehncraft_fv::amg::dist`), the per-rank part of a hierarchy and
its exchange plans, worked out without communication:

```rust
let a = be.amg_matrix()?;                       // gathered, global order
let h = amg::setup::build(a, &settings.amg)?;   // the one-rank hierarchy
let mine = amg::dist::distribute(&h, &owner, &order, rank, size);
// mine.plans[l]: own rows first, then ghosts by owner; peers with send/recv lists
```

`CpuBackend::amg_load` does this itself on a backend made with
`CpuBackend::on_ranks`.

[↑ back to top](#meshes)

---

<a id="flow-ranks"></a>
## 234. Solve a flow on ranks (and check it against one rank)

List rank counts under `[run]` of a steady laminar flow case: after the
run the flow is solved again from the same start on each rank count, the
ranks as threads of this process (recursive coordinate bisection, one
halo layer), and the run fails unless every rank's cells are bitwise
those of one rank (MPI step W5, D93):

```toml
[run]
ranks = [2, 3, 5]
```

```
target/release/foehn flow examples/cases/flow/ranks_channel_poly.toml
```

The log says `solve on N ranks: K outer iterations, bitwise equal to one
rank` for every count; `ranks.vtkhdf` holds the solution of the largest
count with the cell field `rank`. Turbulence, scalars, energy, transient
runs and the coupled solver work too (recipe 244).

Across MPI processes (needs MPI as in recipe 4):

```
cargo build --release -p foehncraft-flow --features mpi --example mpi_flow
for n in 1 2 3 4; do mpirun --oversubscribe -n $n target/release/examples/mpi_flow; done
```

In code (`foehncraft_flow::ranks`), one rank's flow on its part: the
`Target` carries the ranks (`Target { ranks: Some(Ranks::new(comm, &local, &part)), .. }`),
the CPU flow device does the rest, and the one-rank solver runs unchanged:

```rust
let part = partition(&mesh, &geo, size, Method::Rcb)?;
let reference = ranks::solve_one(&case, &mesh, iterations, Some(1))?;
let mine = ranks::solve_part(&case, &mesh, &part, comm.clone(), iterations, Some(1))?;
let c = ranks::check(comm.rank(), &mine, &reference);  // cell by cell, bit for bit
assert!(c.bitwise());
```

[↑ back to top](#meshes)

---

<a id="flow-models-ranks"></a>
## 244. Check a turbulent, heated or transient flow on ranks

`[run] ranks` works with the k-ω SST and the k-ε models (wall functions
included), passive scalars, the energy equation of one fluid, transient
runs (checked with PISO and BDF2) and the coupled solver (MPI step
W5, D93):

```toml
[turbulence]
model = "k-omega-sst"

[energy]
initial = 300.0

[run]
ranks = [2, 3, 5]
```

```
target/release/foehn flow examples/cases/flow/ranks_heated_channel_sst.toml
```

Every rank's velocity, pressure, temperature, scalars, k, ω (or ε) and
ν_t must be bitwise those of one rank, and so must the residual history
(`solve on N ranks: K outer iterations, bitwise equal to one rank`);
`ranks.vtkhdf` holds these fields of the largest rank count. A
transient case is checked over its time steps with the fixed step
`[time] time_step` (`steps`, or `end_time` over the step).

The time of the iterations on ranks against threads only, on a channel
of hexahedra of any size (cells in x, y, z and the iterations; every run
is checked bitwise against the first):

```
cargo run --release -p foehncraft-flow --example ranks_scaling -- 200 100 50 10
```

Not on ranks yet, refused with an error that names them: the ζ-f and
elliptic-blending models, the curvature correction and the QCR, the MLP
limiter and quadratic gradients, gravity, turning zones (MRF, sliding),
solid regions (conjugate heat transfer), a `[[schedule]]`, a level sweep
and a restart. Monitors, probes and balances of a run on ranks come with
the run on ranks itself (MPI step W7, `foehn run --ranks N`).

[↑ back to top](#meshes)
