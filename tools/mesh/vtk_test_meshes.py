"""Writes the VTK test meshes of crates/foehncraft-mesh/tests/data/vtk.

The meshes are the mixed Gmsh test mesh of crates/foehncraft-mesh/tests/data
(hexahedra, wedges, quadrilaterals and triangles as boundary faces, regions
as cell data `gmsh:physical`), converted by meshio, an independent writer of
the VTK formats (only the generator needs it, the repository does not):
    pip install meshio
Run from the repository root:
    python3 -m tools.mesh.vtk_test_meshes
Writes channel_mixed_ascii.vtk, channel_mixed_binary.vtk (legacy),
channel_mixed_ascii.vtu, channel_mixed_binary.vtu (inline base64, no
compression) and channel_mixed_zlib.vtu (zlib).

channel_mixed_appended.vtu (appended raw binary, every cell a polyhedron,
the names as field data) is written by Foehncraft's own writer,
`format::vtk::write_vtu` with `binary = true`, from the same Gmsh mesh.
"""

import meshio

SRC = "crates/foehncraft-mesh/tests/data/channel_mixed_ascii.msh"
OUT = "crates/foehncraft-mesh/tests/data/vtk"


def main():
    m = meshio.read(SRC)
    m.field_data = {}
    m.cell_sets = {}
    m.point_data = {}
    # one block per cell type, the region array as cell data
    m.cell_data = {"gmsh:physical": m.cell_data["gmsh:physical"]}
    meshio.write(f"{OUT}/channel_mixed_ascii.vtk", m, binary=False)
    meshio.write(f"{OUT}/channel_mixed_binary.vtk", m, binary=True)
    meshio.write(f"{OUT}/channel_mixed_ascii.vtu", m, binary=False)
    meshio.write(f"{OUT}/channel_mixed_binary.vtu", m, binary=True, compression=None)
    meshio.write(f"{OUT}/channel_mixed_zlib.vtu", m, binary=True, compression="zlib")


if __name__ == "__main__":
    main()
