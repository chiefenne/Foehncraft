"""Foehncraft from Python (D61): cases, expressions, runs, results, materials
and parametric studies on the same engine as the `foehn` command line.

    from foehncraft import Case, run
    from foehncraft.expr import y

    case = Case.load("examples/cases/flow/poiseuille_hex.toml")
    case["boundary.xmin.velocity"] = [6 * y * (1 - y), 0.0, 0.0]
    del case["boundary.xmin.profile"]
    res = run(case, out="output/poiseuille_py")
    print(res.status, res.orders)

The package works with the native module (built with maturin, runs in
process) or with the `foehn` executable alone (`backend="cli"`).
"""

from . import expr, materials, study
from .case import Case
from .results import Result, read_csv
from .run import FoehncraftError, check, find_foehn, native_available, run
from .tools import call, commands

__all__ = [
    "Case",
    "FoehncraftError",
    "Result",
    "call",
    "check",
    "commands",
    "expr",
    "find_foehn",
    "materials",
    "native_available",
    "read_csv",
    "run",
    "study",
]


def _version() -> str:
    try:
        from ._native import version

        return version()
    except ImportError:
        from importlib.metadata import PackageNotFoundError, version

        try:
            return version("foehncraft")
        except PackageNotFoundError:
            return "unknown"


__version__ = _version()
