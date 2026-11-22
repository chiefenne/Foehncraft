# Meshes

Part of the [Example gallery](../gallery.md): the built-in test mesher in 2D and 3D, partitioning, the Gmsh and SU2 readers and flow on poor Gmsh meshes.

| # | Example | Stage | Feature |
|---|---|---|---|
| 3 | [The four test-mesh types in 3D](#mesh-types-3d) | 1 | built-in mesher, grading, perturbation, quality report |
| 4 | [A polyhedral mesh in four parts](#poly-partition) | 1 | partitioning, halo layout, point location |
| 5 | [The four test-mesh types in 2D](#mesh-types-2d) | 1 | one-cell-thick 2D meshes, boundary renaming |
| 6 | [T-junction from Gmsh](#tjunction) | 1 | Gmsh reader, conversion to polyhedra, renumbering |
| 26 | [Flow on raw Gmsh tetrahedra: the cell-quality remediation](#skewed-tets) | 6 | skewed faces, Gauss pressure gradient, limiter, divergence report (D50, D63) |
| 127 | [A PyAero airfoil mesh from an SU2 file](#su2-airfoil) | 6 | SU2 reader, mesh-format layer, 2D files extruded one layer, markers as boundaries (D187) |
| 129 | [A polyhedral mesh from OpenFOAM and CGNS](#foam-cgns-polyhedra) | 6 | OpenFOAM and CGNS readers, polyhedral cells, one mesh through two formats (D188) |
| 134 | [One polyhedral mesh through five formats](#five-formats) | 6 | industry .msh and VTK readers, VTU writer, the same mesh as OpenFOAM, CGNS, .msh and VTU (D191, D192) |
| 138 | [Taylor–Couette sliding case from a two-zone SU2 file](#su2-zones) | 6 | multi-zone SU2 and CGNS, rotor and stator as cell zones, interface boundaries, `[[sliding]]` (D195) |
| 189 | [Which cells can the curvature correction and the QCR trust? The operator check](#operator-check) | 6 | gradient and Hessian on test fields, `MESH004`–`MESH006`, `LIM008`, `operators.csv` (D235) |

---

<a id="mesh-types-3d"></a>
## 3. The four test-mesh types in 3D

<img src="../figures/gallery/mesh_types_3d.png" alt="Hexahedral, tetrahedral, polyhedral and mixed test meshes of a box, cut open" width="800" height="620">

The built-in test mesher (D34) makes four mesh types of the same box
2 × 1 × 1: hexahedra, tetrahedra (six per hexahedron), polyhedra (the dual
of the tetrahedral mesh, one cell per node) and a mixed mesh of hexahedra,
prisms and tetrahedra, with polyhedra where a hexahedron meets the
triangles of its neighbours. All four are clustered towards the walls
y = 0 and y = 1 (middle-to-wall ratio 4) and have perturbed nodes (15 % of
the local interval), so the geometry is never trivially orthogonal. Every
run writes the quality report (`mesh_report.md` and `.json`): closure of
every cell, total volume against the enclosed volume, non-orthogonality,
skewness, face warp, aspect and volume ratio, inverted faces and negative
volumes. All four meshes have the exact volume 2 and closure below 1e-15.

| | |
|---|---|
| Meshes | hex 576, tet 3 456, poly 819, mixed 1 728 cells |
| Case files | `mesh/hex_box.toml`, `mesh/tet_box.toml`, `mesh/poly_box.toml`, `mesh/mixed_box.toml` |

```
cargo build --release -p foehncraft-cli
for t in hex tet poly mixed; do
  target/release/foehn mesh examples/cases/mesh/${t}_box.toml
done
# output/<case>/mesh_report.md, mesh_report.json, mesh.vtkhdf (open in ParaView)
python3 -m tools.mesh.render_mesh output/poly_box/mesh.vtkhdf poly_box.png --cut
```

[↑ back to top](#meshes)

---

<a id="poly-partition"></a>
## 4. A polyhedral mesh in four parts

<img src="../figures/gallery/poly_partition.png" alt="Polyhedral unit cube in four METIS parts" width="800" height="622">

A polyhedral mesh of the unit cube (2 197 cells, nodes perturbed by 20 %)
split into four parts with METIS (D32), one colour per rank. Each part is
laid out as solver ranks will see it: inner cells, then cells touching the
halo, then the halo cells by source rank, each group in reverse
Cuthill-McKee order; faces inner, processor, then boundary by boundary. The
report lists the cells, halo cells, processor faces and neighbours of each
rank (imbalance 1.001, 1 230 cut faces). The point-location test of D47
throws 5 000 random points and finds each one by walking from a random cell
through the faces; all 5 000 are found in the right cell.

| | |
|---|---|
| Mesh | poly, 12 × 12 × 12 base grid, 2 197 cells |
| Case file | `mesh/poly_partition.toml` |

```
target/release/foehn mesh examples/cases/mesh/poly_partition.toml
target/release/foehn mesh examples/cases/mesh/poly_partition.toml --method rcb --parts 6
python3 -m tools.mesh.render_mesh output/poly_partition/mesh.vtkhdf part.png --cut --field rank --categorical
```

[↑ back to top](#meshes)

---

<a id="mesh-types-2d"></a>
## 5. The four test-mesh types in 2D

<img src="../figures/gallery/mesh_types_2d.png" alt="Hexahedral, tetrahedral, polyhedral and mixed 2D meshes" width="800" height="460">

2D meshes are one cell thick in z (D3): quadrilaterals become hexahedra,
triangles become prisms, polygons become prismatic polyhedra. The mixed
mesh combines both. Clustering towards x = 0 and towards the walls, nodes
perturbed in x and y. The case files rename the generator's boundaries
(`xmin` to `inlet` and so on) with a `[mesh.boundaries]` table; `zmin` and
`zmax` take the symmetry condition in a 2D case.

| | |
|---|---|
| Meshes | hex 160, tet 320, poly 187, mixed 240 cells |
| Case files | `mesh/hex_2d.toml`, `mesh/tet_2d.toml`, `mesh/poly_2d.toml`, `mesh/mixed_2d.toml` |

```
for t in hex tet poly mixed; do
  target/release/foehn mesh examples/cases/mesh/${t}_2d.toml
done
python3 -m tools.mesh.render_mesh output/poly_2d/mesh.vtkhdf poly_2d.png --top
```

[↑ back to top](#meshes)

---

<a id="tjunction"></a>
## 6. T-junction from Gmsh

<img src="../figures/gallery/tjunction_tet.png" alt="T-junction, tetrahedral mesh from Gmsh, cut at the mid plane, four ranks" width="800" height="533">

<img src="../figures/gallery/tjunction_poly.png" alt="The same T-junction converted into polyhedra" width="800" height="533">

The T-junction mixer of stage 7 (D45) meshed with Gmsh
(`examples/meshes/tjunction.geo`, main pipe 100 mm diameter, branch 50 mm)
and read by the Gmsh v4 reader (D33). Physical groups become boundaries
(`inlet_main`, `inlet_branch`, `outlet`, `wall`) and the cell zone `fluid`.
The second case converts the 77 776 tetrahedra into 15 927 polyhedra, one
per node; on the curved wall the boundary faces follow the triangles of
the tetrahedral mesh, so both meshes have the same volume to the last
digit (8.310580753620e-3 m³). Both are renumbered by reverse
Cuthill-McKee and split into four parts. The point-location walk stops at
the curved, non-convex wall for some points (176 of 1 195); the bin search
finds those, none lands in a wrong cell.

| | |
|---|---|
| Meshes | tet 77 776 cells, poly 15 927 cells (7 to 23 faces, mean 15) |
| Case files | `mesh/tjunction_tet.toml`, `mesh/tjunction_poly.toml` |

```
gmsh -3 examples/meshes/tjunction.geo -format msh41 -o examples/meshes/tjunction.msh
target/release/foehn mesh examples/cases/mesh/tjunction_tet.toml
target/release/foehn mesh examples/cases/mesh/tjunction_poly.toml
python3 -m tools.mesh.render_mesh output/tjunction_poly/mesh.vtkhdf tj.png \
    --field rank --categorical --half --top --crop 0.2 0.8 -0.06 0.31
```

[↑ back to top](#meshes)

---

<a id="skewed-tets"></a>
## 26. Flow on raw Gmsh tetrahedra: the cell-quality remediation

<img src="../figures/gallery/pipe_gmsh_raw_residuals.svg" alt="Residuals of the pipe on raw Gmsh tetrahedra with and without the remediation" width="800" height="364">

Gmsh's tetrahedra without its optimisation pass contain faces of a
skewness above 0.9: the face centroid lies almost as far from the line
through the two cell centroids as the centroids lie apart, mostly in
tetrahedra with two faces on a curved wall. On such faces the
least-squares pressure gradient (in the momentum source, the Rhie–Chow
flux and the wall pressure) drives the steady SIMPLEC iteration away
from the solution. The pipe of `examples/meshes/pipe_gmsh_raw.geo`
(D = 0.1 m, L = 1 m, Re = 100, 1 572 tetrahedra, 28 faces above 0.9, the
largest 1.23) shows it: with `solver.remediation = "off"` the residuals
stall at 0.1, and in iteration 91 the pressure-correction matrix stops
being positive definite. The divergence report (D63) names the skewed
faces as the first cause (`DIV003`) with the fix
`solver.remediation="auto"`.

With the default `solver.remediation = "auto"` the solver measures the
faces at start-up, reports them (`MESH001`), switches to the Gauss
pressure gradient (the face sum of the face pressures, with the boundary
pressures) and the Venkatakrishnan limiter on the velocity gradients, says
so (`MESH002`, `flow.log`, `flow.md`) and converges in 148 iterations (105 before the bounded face velocity of
2026-11-20, the run of the chart):

```
mesh quality: 28 faces with skewness above 0.9 (largest 1.233) in 44 cells; an optimised mesh avoids them (Gmsh: Mesh.OptimizeNetgen = 1)
cell-quality remediation (D50, solver.remediation = "auto"): solver.pressure_gradient: "least-squares" -> "gauss", solver.momentum.limiter: "none" -> "venkatakrishnan", solver.face_velocity: "interpolated" -> "bounded"
```

**Bounded face velocity (2026-11-20).** Turbulent MRF runs on raw Gmsh
tetrahedra (the stirred tank of the Engineering Gallery) still diverged
with the Gauss gradient and the limiter: the velocity in a few cells at
the blade–hub corner grew about twofold every 50 iterations until the
torque ran away. The growing term is the face velocity of the Rhie–Chow
flux, w u_P + (1 − w) u_N + ∇ū_f·(x_f − x_w): on a skewed face the shift
x_f − x_w is as long as the centroid distance, and the unlimited
gradients feed a velocity mode back into continuity with a gain above
one. The remediation now takes the mean of the two limited
extrapolations, w (u_P + ψ_P ∇u_P·r_P) + (1 − w)(u_N + ψ_N ∇u_N·r_N),
which is just as exact for a linear field and stays between the values
the two cells see (`solver.face_velocity = "bounded"`). On the ungraded
58k-cell tank this removes the blow-up (stable past iteration 500, torque
−1.03 N·m; before: runaway from iteration 330). A coarser version of the
tank (about 36k tetrahedra, faces up to skewness 1.3) is the test case:
`python3 examples/studies/skewed_tets.py` meshes it and runs both face
velocities for 300 iterations, a few minutes each:

```
iteration | bounded: residual     torque | interpolated: residual     torque
      100 |          1.20e-02      -1.13 |               1.18e-02      -1.09
      200 |          4.08e-03      -1.09 |               7.71e-02      -1.13
      300 |          1.73e-03      -1.11 |               1.23e-01       -177
```

The switch is global, not per cell: treating only the cells at the
skewed faces (first-order convection, a Gauss gradient or a damped
correction there) did not converge on these meshes, because the error of
the pressure gradient travels with the pressure correction through the
whole domain. The Gauss pressure gradient keeps second order on
tetrahedra (observed 2.00 on the Poiseuille sweep), but it is less
accurate than least squares on smooth meshes (the lift of the cylinder
at Re = 20 is 3.5 % lower), so it stays a remedy for poor meshes and
least squares the default. A mesh without faces above
`solver.skewness_limit` (0.9) runs bitwise as before; none of the
catalogue cases has such faces.

The better cure is the mesh: `Mesh.OptimizeNetgen = 1` in the .geo file
(as in `examples/meshes/tjunction_netgen.geo`) removes the skewed faces.
The T-junction of gallery 6 has 461 faces above 0.9 in Gmsh's raw mesh
and none after the Netgen pass; a laminar run on the raw mesh diverged
without the remediation and converges with it, and the optimised mesh
runs unchanged.

| | |
|---|---|
| Mesh | `examples/meshes/pipe_gmsh_raw.msh`, 1 572 tetrahedra (in the repository) |
| Case file | `flow/pipe_gmsh_raw.toml` |
| Keys | `solver.remediation`, `solver.skewness_limit`, `solver.pressure_gradient`, `solver.face_velocity`, `solver.transposed_stress` |

```
cargo build --release -p foehncraft-cli
target/release/foehn flow examples/cases/flow/pipe_gmsh_raw.toml --set solver.log_every=1
target/release/foehn flow examples/cases/flow/pipe_gmsh_raw.toml --set solver.remediation=off \
    --set solver.log_every=1 --out output/pipe_gmsh_raw_off
python3 -m tools.plot.residuals_chart pipe_gmsh_raw_residuals.svg \
    "output/pipe_gmsh_raw_off=remediation off" "output/pipe_gmsh_raw=remediation auto"
```

[↑ back to top](#meshes)

---

<a id="su2-airfoil"></a>
## 127. A PyAero airfoil mesh from an SU2 file

<img src="../figures/gallery/naca0012_su2_speed.png" alt="Speed around a NACA 0012 at Re 1000 and 4 degrees on the PyAero C-mesh read from an SU2 file" width="800" height="400">

<img src="../figures/gallery/naca0012_su2_vorticity.png" alt="Vorticity of the same flow: the boundary layers and the wake" width="800" height="400">

<img src="../figures/gallery/naca0012_su2_mesh.png" alt="The PyAero C-mesh around the airfoil, a blunt trailing edge with its own block" width="800" height="400">

A NACA 0012 C-mesh made by [PyAero](https://github.com/chiefenne/PyAero)
in batch mode (27 297 quadrilaterals, orthogonal near the wall, a blunt
trailing edge of 0.4 % chord with its own block) and written in the SU2
format. Foehncraft reads it through the mesh-format layer (D187): the file
is 2D (`NDIME= 2`), so the reader extrudes it one cell layer in z, the
quadrilaterals become hexahedra and the planes the boundaries `zmin` and
`zmax`; the five markers `airfoil`, `inlet`, `outlet`, `top` and `bottom`
become boundaries by their names. The quality report finds the mesh valid
(every cell closed to 1e-16, the volume equal to the volume the boundary
encloses to 4e-16, non-orthogonality at most 35°, skewness at most 0.14).
The flow is laminar at Re = 1000 and 4° angle of attack (U = 1 m/s,
chord 1 m, ν = 1e-3 m²/s); the coupled solver converges in 56 iterations
and 10 s on four cores, with C_D = 0.128 and C_L = 0.209 (along and across
the free stream). The repository holds a coarser PyAero mesh of the same
airfoil (6 420 cells, `examples/meshes/naca0012.su2`) that the case file
uses and CI runs in about a second: C_D = 0.128, C_L = 0.213. The far field
is only 3.5 chords above and below the airfoil, so the run shows the
reader, not a validated lift.

| | |
|---|---|
| Meshes | PyAero: 27 297 hexahedra after the extrusion (pictures), 6 420 in the repository (case, CI) |
| Case files | `mesh/naca0012_su2.toml`, `flow/naca0012_su2_re1000.toml` |
| PyAero control files | `examples/meshes/naca0012_pyaero_fine.json` (pictures), `naca0012_pyaero.json` (the mesh in the repository) |

```
# in a PyAero clone (https://github.com/chiefenne/PyAero):
python src/PyAero.py -no-gui /path/to/Foehncraft/examples/meshes/naca0012_pyaero_fine.json
cp data/OUTPUT/NACA0012.su2 /path/to/Foehncraft/examples/meshes/naca0012_fine.su2
# in Foehncraft:
target/release/foehn mesh examples/meshes/naca0012_fine.su2
target/release/foehn flow examples/cases/flow/naca0012_su2_re1000.toml \
    --set mesh.file=../../meshes/naca0012_fine.su2 --out output/naca0012_fine
python3 -m tools.mesh.render_mesh output/naca0012_fine/fields.vtkhdf speed.png \
    --field speed --top --crop -1.2 3.5 -1.5 1.5 --no-edges --size 1200x600 --zoom 2.3
python3 -m tools.mesh.render_mesh output/naca0012_fine/fields.vtkhdf vorticity.png \
    --field vorticity_z --diverging --range -20 20 --top --crop -1.2 3.5 -1.5 1.5 \
    --no-edges --size 1200x600 --zoom 2.3
python3 -m tools.mesh.render_mesh output/naca0012_fine/fields.vtkhdf mesh.png \
    --top --crop -0.4 1.6 -0.6 0.6 --size 1200x600 --zoom 2.0
```

<a id="foam-cgns-polyhedra"></a>
## 129. A polyhedral mesh from OpenFOAM and CGNS

The dual of a tetrahedral unit cube (OpenFOAM's `polyDualMesh`, 81
polyhedral cells with 9 to 48 faces each) is read from the OpenFOAM
`polyMesh` directory and from a CGNS file of the same mesh (NGON_n faces and
NFACE_n cells, HDF5). Both readers give the same 81 cells and the same six
boundaries `xmin` to `zmax` with 40 faces each; the quality report finds the
mesh valid (volume 1 and the volume enclosed by the boundary equal to
4e-16, every boundary of area 1). The same two readers take hexahedral,
tetrahedral and prism meshes: the repository holds Gmsh CGNS files of a unit
cube in each type, and an OpenFOAM channel in ASCII and in binary that reads
alike.

```
target/release/foehn mesh crates/foehncraft-mesh/tests/data/openfoam/dual_polyhedra
target/release/foehn mesh crates/foehncraft-mesh/tests/data/cgns/dual_ngon.cgns
```

| | |
|---|---|
| Meshes | `crates/foehncraft-mesh/tests/data/openfoam/{channel_ascii,channel_binary,dual_polyhedra}`, `.../cgns/{box_hex,box_tet,box_prism,dual_ngon}.cgns` |
| Made by | OpenFOAM utilities and Gmsh; `tools/mesh/cgns_test_meshes.py`, `tools/mesh/cgns_ngon_fixture.py`, `tools/mesh/openfoam_test_meshes.md` |
| Recipe | [recipe 129](../recipes/meshes.md#mesh-cgns-openfoam) |

<a id="five-formats"></a>
## 134. One polyhedral mesh through five formats

The dual of a tetrahedral unit cube (81 polyhedral cells, 582 faces, made by
OpenFOAM's `polyDualMesh`) exists in the repository as an OpenFOAM
`polyMesh` (gallery 129), as a CGNS NGON/NFACE file, as an industry `.msh`
(`tests/data/industry_msh/dual_polyhedra.msh`) and, written by Foehncraft
itself, as a VTK XML grid. Every reader gives the same 81 cells, the same
six boundaries `xmin` to `zmax` of area 1 and the volume 1 (4e-16 closed);
a `.vtu` written with `output.vtu` reads back with the same names. A mixed
Gmsh mesh of hexahedra and wedges converted by an independent tool to
legacy ASCII, legacy binary, XML ASCII, XML base64 and zlib-compressed XML
reads as the mesh of the Gmsh file in all five encodings.

```
target/release/foehn mesh crates/foehncraft-mesh/tests/data/industry_msh/dual_polyhedra.msh --vtu
target/release/foehn mesh output/dual_polyhedra/mesh.vtu
```

| | |
|---|---|
| Meshes | `crates/foehncraft-mesh/tests/data/industry_msh/*.msh`, `.../vtk/channel_mixed_*.{vtk,vtu}` |
| Made by | `tools/mesh/industry_msh_fixtures.py` (from the OpenFOAM fixtures), `tools/mesh/vtk_test_meshes.py` (meshio, generator only) |
| Recipe | [recipe 134](../recipes/meshes.md#mesh-msh-vtk) |

[↑ back to top](#meshes)

<a id="su2-zones"></a>
## 138. Taylor–Couette sliding case from a two-zone SU2 file

The gap between two coaxial cylinders, written the way a mesher writes a
machine: an SU2 file with two zones, the rotor (radius 0.05 to 0.075 m, 32
sectors) and the stator (0.075 to 0.10 m, 48 sectors), no shared nodes and an
interface that does not match. The reader puts both zones in one mesh (cell
zones `zone1` and `zone2`), the markers `rotor_interface` and
`stator_interface` are the two boundaries, and the sliding interface imprints
them at the rotor angle in every step, as in the generated case of the
imprinting study (same flow, Re = 50, same check against the analytic
Taylor–Couette solution). The run takes 4 s on 640 cells: the relative L2
error of the velocity is 1.3e-2 against the band of 3e-2.

```
target/release/foehn flow examples/cases/flow/taylor_couette_su2_zones.toml
```

| | |
|---|---|
| Case | `examples/cases/flow/taylor_couette_su2_zones.toml` |
| Mesh | `examples/meshes/taylor_couette_zones.su2`, made by `tools/mesh/taylor_couette_zones.py` |
| Recipe | [recipe 138](../recipes/meshes.md#mesh-multizone) |

[↑ back to top](#meshes)

---

<a id="operator-check"></a>
## 189. Which cells can the curvature correction and the QCR trust? The operator check

<img src="../figures/gallery/operator_check.png" alt="Three panels of a shallow 2:1 cavity, cells coloured on a log scale from 0.001 (black) to 2 (light yellow). Top: regular hexahedra, the Hessian error is zero inside and about 0.5 in the two layers of cells along the walls (368 of 2048 cells flagged). Middle: perturbed tetrahedra, the gradient error is about 0.2 to 0.4 in every cell. Bottom: the same tetrahedra, the Hessian error is about 0.5 to 1 in every cell" width="900" height="1150">

**What is checked.** Before the solve the run applies its velocity
gradient (least squares over the face neighbours, D24) and the gradient
of that gradient (the velocity Hessian behind DS_ij/Dt of the curvature
correction) to linear and quadratic test fields on its own mesh. The
gradient must be exact for a linear field on every cell; where it is also
exact for a quadratic field it is second order, and only there do the
operators that differentiate it once more (the divergence of the explicit
QCR stress, DS_ij/Dt) converge with the mesh. The errors per cell are
dimensionless: the gradient error over h |H| (h the longest stencil
distance), the Hessian error over |H|.

**Results.** On regular hexahedra (top) both operators are exact inside;
the cells at the walls have one-sided stencils, so the gradient is first
order in one layer and the Hessian in two. On tetrahedra with perturbed
vertices (middle, bottom) every cell is flagged: gradient errors of 0.2 to
0.4, Hessian errors near one. This is the reason why the curvature
correction and the QCR did not converge on tetrahedra in the manufactured
solutions of #184 (gallery 184). The run reports it before it starts:
`MESH005` (curvature correction) and `MESH006` (QCR) with the count and the
worst cells, `LIM008` in `result.json`, the 50 worst cells per operator in
`flow.json` and `foehn flow --check`, and the errors of every cell in
`operators.csv`. Nothing is reported where every operator the run uses is
exact, and the check changes no result. The next steps of D235 treat the
flagged cells (a quadratic least-squares fit over the neighbours of the
neighbours for the QCR gradient and the Hessian, then a local remesh of
the cells still flagged).

| | |
|---|---|
| Case file | `flow/cavity_sst_operator_check.toml` |
| Plot | `tools/plot/operator_check.py` |
| Recipe | [189](../recipes/meshes.md#operator-check) |

```
F=examples/cases/flow/cavity_sst_operator_check.toml
target/release/foehn flow $F --out output/ops_tet --set output.vtkhdf=true
target/release/foehn flow $F --out output/ops_hex --set output.vtkhdf=true \
    --set mesh.generate.type=\"hex\" --set mesh.generate.perturb=0
python3 -m tools.plot.operator_check docs/figures/gallery/operator_check.png output/ops_hex output/ops_tet
```

[↑ back to top](#meshes)
