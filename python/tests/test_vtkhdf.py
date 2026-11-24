"""`foehncraft.vtkhdf.internal`: the blocks of a VTKHDF file in one piece."""

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")

from foehncraft.vtkhdf import internal  # noqa: E402


def write_blocked(path):
    """Two tetrahedra A = (0, 1, 2, 3) and B = (1, 2, 3, 4) sharing a face,
    each owned by its own block with the other as a ghost; block 0 owns B,
    and its points come in reverse order of their identities."""
    xyz = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]], dtype=float)
    with h5py.File(path, "w") as f:
        g = f.create_group("VTKHDF/internal")
        ids0 = np.array([4, 3, 2, 1, 0])  # block 0: B then the ghost A
        ids1 = np.array([0, 1, 2, 3, 4])  # block 1: A then the ghost B
        local0 = {p: k for k, p in enumerate(ids0)}
        local1 = {p: k for k, p in enumerate(ids1)}
        a, b = [0, 1, 2, 3], [1, 2, 3, 4]
        conn = [local0[p] for p in b + a] + [local1[p] for p in a + b]
        g["NumberOfCells"] = [2, 2]
        g["NumberOfPoints"] = [5, 5]
        g["NumberOfConnectivityIds"] = [8, 8]
        g["Points"] = np.concatenate([xyz[ids0], xyz[ids1]])
        g["Types"] = np.full(4, 10, dtype=np.uint8)
        g["Offsets"] = [0, 4, 8, 0, 4, 8]
        g["Connectivity"] = conn
        cd = g.create_group("CellData")
        cd["pressure"] = [2.0, 1.0, 1.0, 2.0]
        cd["vtkGhostType"] = np.array([0, 1, 0, 1], dtype=np.uint8)
        cd["GlobalCellIds"] = [0, 1, 1, 0]
        cd["OriginalCellIds"] = [1, 0, 0, 1]
        pd = g.create_group("PointData")
        pd["GlobalPointIds"] = np.concatenate([ids0, ids1])
        pd["vtkGhostType"] = np.array([0, 0, 0, 0, 1, 0, 1, 1, 1, 1], dtype=np.uint8)
    return xyz


def test_blocks_come_back_in_one_piece(tmp_path):
    p = tmp_path / "blocked.vtkhdf"
    xyz = write_blocked(p)
    g = internal(p)
    np.testing.assert_array_equal(g["Points"], xyz)
    np.testing.assert_array_equal(g["Types"], [10, 10])
    np.testing.assert_array_equal(g["Offsets"], [0, 4, 8])
    np.testing.assert_array_equal(g["Connectivity"], [0, 1, 2, 3, 1, 2, 3, 4])
    np.testing.assert_array_equal(g["CellData/pressure"], [1.0, 2.0])
    assert sorted(g["CellData"]) == ["pressure"]
    assert g["PointData"] == {}


def test_a_file_without_blocks_comes_back_as_it_is(tmp_path):
    p = tmp_path / "whole.vtkhdf"
    with h5py.File(p, "w") as f:
        g = f.create_group("VTKHDF/internal")
        g["NumberOfCells"] = [1]
        g["NumberOfPoints"] = [4]
        g["NumberOfConnectivityIds"] = [4]
        g["Points"] = np.eye(4, 3)
        g["Types"] = np.array([10], dtype=np.uint8)
        g["Offsets"] = [0, 4]
        g["Connectivity"] = [3, 2, 1, 0]
        g.create_group("CellData")["zone"] = [0.0]
    with h5py.File(p, "r") as f:
        g = internal(f)
    np.testing.assert_array_equal(g["Points"], np.eye(4, 3))
    np.testing.assert_array_equal(g["Connectivity"], [3, 2, 1, 0])
    np.testing.assert_array_equal(g["CellData"]["zone"], [0.0])
