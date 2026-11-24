"""Reads Foehncraft VTKHDF files with VTK's own reader (the reader ParaView
uses) and checks them:

- the file loads, with the block `internal` and one block per boundary;
- the cell volumes agree with the field `volume` written by Foehncraft. For
  standard cells VTK's own cell volumes are used (this checks the node order
  of tetrahedra, hexahedra, prisms and pyramids); for polyhedra the volume is
  computed here from the faces as VTK reads them, by the divergence theorem
  (this checks that every face is there and points outward; VTK's own
  polyhedron volume is the volume of the convex hull of the nodes, which
  differs for warped faces).

    pip install vtk
    python3 -m tools.check.vtkhdf FILE.vtkhdf [...]
"""
import sys
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy as v2n


def polyhedron_volume(cell):
    vol = 0.0
    for i in range(cell.GetNumberOfFaces()):
        f = cell.GetFace(i)
        p = np.array([f.GetPoints().GetPoint(k) for k in range(f.GetNumberOfPoints())])
        m = p.mean(axis=0)
        for k in range(len(p)):
            q = p[(k + 1) % len(p)]
            s = 0.5 * np.cross(p[k] - m, q - m)
            vol += np.dot((p[k] + q + m) / 3.0, s) / 3.0
    return vol


def check(path):
    r = vtk.vtkHDFReader()
    r.SetFileName(path)
    r.Update()
    out = r.GetOutputDataObject(0)
    assert out.IsA("vtkPartitionedDataSetCollection"), out.GetClassName()
    names = [out.GetMetaData(i).Get(vtk.vtkCompositeDataSet.NAME())
             for i in range(out.GetNumberOfPartitionedDataSets())]
    parts = out.GetPartitionedDataSet(names.index("internal"))
    n, types, err = 0, set(), 0.0
    for k in range(parts.GetNumberOfPartitions()):
        grid = parts.GetPartition(k)
        # a blocked file (layout 1) repeats the cells next to a block as
        # ghosts; each cell is checked in the block that owns it
        ghost = grid.GetCellData().GetArray("vtkGhostType")
        own = np.arange(grid.GetNumberOfCells())
        if ghost is not None:
            own = np.nonzero(v2n(ghost) == 0)[0]
        n += len(own)
        t = np.array([grid.GetCellType(int(i)) for i in own], dtype=int)
        types |= set(t.tolist())
        ours = grid.GetCellData().GetArray("volume")
        # Volumes are compared only where Foehncraft wrote them (the mesh
        # files); computing them for a field file without `volume` checked
        # nothing and took minutes for polyhedral meshes.
        if ours is None or len(own) == 0:
            continue
        sizes = vtk.vtkCellSizeFilter()
        sizes.SetInputData(grid)
        sizes.Update()
        vol = v2n(sizes.GetOutput().GetCellData().GetArray("Volume"))[own].copy()
        for j in np.nonzero(t == vtk.VTK_POLYHEDRON)[0]:
            vol[j] = polyhedron_volume(grid.GetCell(int(own[j])))
        ours = v2n(ours)[own]
        err = max(err, abs(vol - ours).max() / abs(ours).max())
    msg = f"{path}: {n} cells in {parts.GetNumberOfPartitions()} block(s), {len(names) - 1} boundaries"
    if err > 0.0 or any(
        parts.GetPartition(k).GetCellData().GetArray("volume") is not None
        for k in range(parts.GetNumberOfPartitions())
    ):
        msg += f", max relative volume difference {err:.1e}"
        assert err < 1e-8, msg
    msg += f", VTK cell types {sorted(types)}"
    print(msg)


for p in sys.argv[1:]:
    check(p)
