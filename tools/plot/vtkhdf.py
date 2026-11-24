"""The block `internal` of a FoehnCraft VTKHDF file in one piece, for the
plot tools: `fields.vtkhdf` is written in blocks with ghost cells (blocked
layout 1), and `internal(f)` returns the cells in the order of their
identities, each once, as the arrays of one unstructured grid. It is
`foehncraft.vtkhdf.internal` of the repository's Python package (h5py and
NumPy only; the native module is not needed).
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "python"))

from foehncraft.vtkhdf import internal  # noqa: E402

__all__ = ["internal"]
