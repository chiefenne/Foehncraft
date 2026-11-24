"""The cells of a FoehnCraft VTKHDF file in one piece (needs h5py and numpy).

`fields.vtkhdf` is written in blocks (the blocked layout 1 that MurmelCraft
and KoflCraft share): the cells in partitions along a Hilbert curve, each
partition with its own points and a layer of ghost cells from the next
blocks. ParaView and MurmelCraft read that as it is. A script that reads the
arrays itself gets them back in one piece from `internal`:

    import h5py
    from foehncraft.vtkhdf import internal

    with h5py.File("output/run/fields.vtkhdf", "r") as f:
        g = internal(f)
    g["Points"], g["Connectivity"], g["Offsets"], g["Types"]
    g["CellData"]["pressure"]          # or g["CellData/pressure"]

The cells come in the order of their identities, each once (ghosts
dropped), with the cell data in the same order; the points come in the order
of their identities, each once. Polyhedra keep `FaceOffsets`,
`FaceConnectivity`, `PolyhedronOffsets` and `PolyhedronToFaces`. A file of
one partition without ghosts (written before the blocks) comes back as it
is.
"""

from __future__ import annotations

# the arrays of the layout itself, which a script of the whole grid does not want
LAYOUT_CELL = ("vtkGhostType", "GlobalCellIds", "OriginalCellIds")
LAYOUT_POINT = ("vtkGhostType", "GlobalPointIds")


class Grid(dict):
    """Arrays by name, as an h5py group: `g["CellData/p"]` is `g["CellData"]["p"]`."""

    def __getitem__(self, key):
        if isinstance(key, str) and "/" in key:
            head, rest = key.split("/", 1)
            return dict.__getitem__(self, head)[rest]
        return dict.__getitem__(self, key)


def _gather(data, starts, lens):
    """The runs `data[s:s + n]` for each start and length, one after the other."""
    import numpy as np

    lens = np.asarray(lens, dtype=np.int64)
    total = int(lens.sum())
    if total == 0:
        return data[:0]
    ends = np.cumsum(lens)
    idx = np.arange(total, dtype=np.int64) - np.repeat(ends - lens, lens)
    return data[np.repeat(np.asarray(starts, dtype=np.int64), lens) + idx]


def _runs(offsets, counts, base=None):
    """Start and length of each row's run in a concatenated array, from the
    per-partition `offsets` (`counts[p] + 1` values each)."""
    import numpy as np

    starts, lens, row = [], [], 0
    base = base if base is not None else np.zeros(len(counts), dtype=np.int64)
    for p, n in enumerate(counts):
        o = offsets[row : row + n + 1]
        starts.append(base[p] + o[:-1])
        lens.append(np.diff(o))
        row += n + 1
    if not starts:
        return np.zeros(0, np.int64), np.zeros(0, np.int64)
    return np.concatenate(starts), np.concatenate(lens)


def internal(f) -> Grid:
    """The block `internal` of an open h5py file (or of the file at path `f`)
    as one grid; see the module documentation."""
    import h5py
    import numpy as np

    if not isinstance(f, (h5py.File, h5py.Group)):
        with h5py.File(f, "r") as h:
            return internal(h)
    g = f["VTKHDF/internal"] if "VTKHDF" in f else f
    arr = lambda name, d=None: np.asarray(g[name][()]) if name in g else d  # noqa: E731
    nc = arr("NumberOfCells").astype(np.int64).ravel()
    npt = arr("NumberOfPoints").astype(np.int64).ravel()
    nk = arr("NumberOfConnectivityIds").astype(np.int64).ravel()
    cell = {k: np.asarray(v[()]) for k, v in g["CellData"].items()} if "CellData" in g else {}
    point = {k: np.asarray(v[()]) for k, v in g["PointData"].items()} if "PointData" in g else {}
    pts = arr("Points")
    pstart = np.concatenate([[0], np.cumsum(npt)[:-1]]).astype(np.int64)
    kstart = np.concatenate([[0], np.cumsum(nk)[:-1]]).astype(np.int64)

    # points: each once, in the order of their identities
    gpid = point.get("GlobalPointIds")
    if gpid is not None:
        ids, first = np.unique(gpid, return_index=True)
        points = pts[first]
        to_global = np.searchsorted(ids, gpid)
    else:
        points = pts
        to_global = np.arange(len(pts), dtype=np.int64)

    def global_ids(flat, counts):
        """Partition-local point ids of a concatenated array as global ones."""
        part = np.repeat(np.arange(len(counts)), counts)
        return to_global[pstart[part] + flat]

    conn = global_ids(arr("Connectivity").astype(np.int64), nk)
    cs, cl = _runs(arr("Offsets").astype(np.int64), nc, kstart)

    # cells: the owned ones, in the order of their identities
    keep = np.ones(int(nc.sum()), dtype=bool)
    if "vtkGhostType" in cell:
        keep = cell["vtkGhostType"] == 0
    order = np.nonzero(keep)[0]
    if "OriginalCellIds" in cell:
        order = order[np.argsort(cell["OriginalCellIds"][order], kind="stable")]
    types = arr("Types")[order]
    lens = cl[order]
    out = Grid()
    out["NumberOfCells"] = np.array([len(order)], dtype=np.int64)
    out["NumberOfPoints"] = np.array([len(points)], dtype=np.int64)
    out["Points"] = points
    out["Types"] = types
    out["Offsets"] = np.concatenate([[0], np.cumsum(lens)]).astype(np.int64)
    out["Connectivity"] = _gather(conn, cs[order], lens)
    out["NumberOfConnectivityIds"] = np.array([len(out["Connectivity"])], dtype=np.int64)

    if "PolyhedronOffsets" in g:
        nf = arr("NumberOfFaces").astype(np.int64).ravel()
        nfk = arr("NumberOfFaceConnectivityIds").astype(np.int64).ravel()
        fstart = np.concatenate([[0], np.cumsum(nf)[:-1]]).astype(np.int64)
        fkstart = np.concatenate([[0], np.cumsum(nfk)[:-1]]).astype(np.int64)
        fconn = global_ids(arr("FaceConnectivity").astype(np.int64), nfk)
        fs, fl = _runs(arr("FaceOffsets").astype(np.int64), nf, fkstart)
        # the faces of each cell, as rows of all faces
        ps, pl = _runs(arr("PolyhedronOffsets").astype(np.int64), nc)
        part = np.repeat(np.arange(len(nc)), nc)
        nfaces = int(nf.sum())
        to_faces = arr("PolyhedronToFaces", np.zeros(0, np.int64)).astype(np.int64)
        ntf = arr("NumberOfPolyhedronToFaceIds", nf).astype(np.int64).ravel()
        tstart = np.concatenate([[0], np.cumsum(ntf)[:-1]]).astype(np.int64)
        faces = _gather(to_faces, tstart[part[order]] + ps[order], pl[order])
        faces = faces + np.repeat(fstart[part[order]], pl[order])
        flen = fl[faces] if nfaces else np.zeros(0, np.int64)
        out["NumberOfFaces"] = np.array([len(faces)], dtype=np.int64)
        out["PolyhedronOffsets"] = np.concatenate([[0], np.cumsum(pl[order])]).astype(np.int64)
        out["PolyhedronToFaces"] = np.arange(len(faces), dtype=np.int64)
        out["FaceOffsets"] = np.concatenate([[0], np.cumsum(flen)]).astype(np.int64)
        out["FaceConnectivity"] = _gather(fconn, fs[faces] if nfaces else fs[:0], flen)
        out["NumberOfFaceConnectivityIds"] = np.array(
            [len(out["FaceConnectivity"])], dtype=np.int64
        )
        out["NumberOfPolyhedronToFaceIds"] = np.array([len(faces)], dtype=np.int64)

    out["CellData"] = {k: v[order] for k, v in cell.items() if k not in LAYOUT_CELL}
    if gpid is not None:
        out["PointData"] = {k: v[first] for k, v in point.items() if k not in LAYOUT_POINT}
    else:
        out["PointData"] = {k: v for k, v in point.items() if k not in LAYOUT_POINT}
    return out
