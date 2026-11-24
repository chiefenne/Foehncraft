"""Results of a run (D61): the summary, the CSV tables and the fields.

    from foehncraft import Result

    res = Result.load("output/cavity_re100")
    res.ok, res.exit_code, res.status
    res.solves[0]["errors"]["l2"]
    hist = res.table("residuals")          # columns as lists
    hist["iteration"], hist["continuity"]
    f = res.fields()                       # needs h5py and numpy
    f["velocity_x"], f["pressure"]

Everything is read from the files of the output directory, so a `Result`
of an earlier run, of a run on a cluster or of the command line is the
same as that of a run in this process.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any


class Result:
    """The output directory of a run and its result summary
    (`result.json`, `foehncraft.result`)."""

    def __init__(self, out, summary: dict | None = None):
        self.out = Path(out)
        if summary is None:
            with open(self.out / "result.json", encoding="utf-8") as f:
                summary = json.load(f)
        self.summary: dict[str, Any] = summary

    @classmethod
    def load(cls, out) -> Result:
        """Reads `out/result.json`."""
        return cls(out)

    # -- the summary
    @property
    def status(self) -> str:
        """`finished`, `failed-checks`, `stopped` or `error`."""
        return self.summary["status"]

    @property
    def exit_code(self) -> int:
        """The exit code of D64: 0, 1 (error), 2 (failed checks), 3 (stopped)."""
        return self.summary["exit_code"]

    @property
    def ok(self) -> bool:
        """Finished with every check passed."""
        return self.exit_code == 0

    @property
    def limitations(self) -> list[dict]:
        """What to know before trusting the result (D199): `code` (LIM001 to
        LIM006), `feature` and `message`. A run with limitations is still
        `finished` with exit code 0."""
        return self.summary.get("limitations", [])

    @property
    def solves(self) -> list[dict]:
        """One entry per level or time step of a sweep."""
        return self.summary.get("solves", [])

    @property
    def orders(self) -> list[float]:
        """Observed orders in space."""
        return self.summary.get("orders", [])

    @property
    def time_orders(self) -> list[float]:
        """Observed orders in time."""
        return self.summary.get("time_orders", [])

    @property
    def scalar_orders(self) -> dict[str, list[float]]:
        """Observed orders of the passive scalars with an exact solution,
        by name (`<name> (time)` for a time-step sweep)."""
        return {name: orders for name, orders in self.summary.get("scalar_orders", [])}

    @property
    def uncertainty(self) -> dict[str, dict]:
        """Discretisation uncertainty of a sweep of three or more meshes or
        time steps (D67), by quantity: `values` (finest first), `order`,
        `extrapolated`, `band` (φ1 ± band), `gci`, `verdict`."""
        return {e["quantity"]: e for e in self.summary.get("uncertainty", [])}

    @property
    def diagnostics(self) -> list[dict]:
        """Errors, failed checks, warnings and notes with their codes."""
        return self.summary.get("diagnostics", [])

    @property
    def files(self) -> dict[str, str]:
        """The output files with their SHA-256 (D84)."""
        return {f["name"]: f["sha256"] for f in self.summary.get("files", [])}

    def value(self, key: str, solve: int = -1, default=None):
        """A value of a solve by dotted key, e.g. `"errors.l2"`,
        `"benchmark.cd_max"`, `"boundaries.outlet.mass_flow"` (list entries
        by their `name` or index)."""
        d: Any = self.solves[solve] if self.solves else {}
        for k in key.split("."):
            if isinstance(d, list):
                named = [x for x in d if isinstance(x, dict) and x.get("name") == k]
                if named:
                    d = named[0]
                    continue
                try:
                    d = d[int(k)]
                except (ValueError, IndexError):
                    return default
            elif isinstance(d, dict) and k in d:
                d = d[k]
            else:
                return default
        return d

    def flat(self, solve: int = -1) -> dict[str, Any]:
        """The numbers of a solve with dotted keys (for study tables)."""
        out: dict[str, Any] = {}

        def walk(d, path):
            if isinstance(d, dict):
                for k, v in d.items():
                    walk(v, f"{path}.{k}" if path else k)
            elif isinstance(d, list):
                for i, v in enumerate(d):
                    # entries with a name (boundaries) by their name
                    key = v["name"] if isinstance(v, dict) and "name" in v else str(i)
                    walk(v, f"{path}.{key}")
            elif isinstance(d, (int, float)) and not isinstance(d, bool):
                out[path] = d

        if self.solves:
            walk(self.solves[solve], "")
        return out

    def __repr__(self) -> str:
        return f"Result({str(self.out)!r}, status={self.status!r})"

    # -- tables
    def table(self, name: str) -> dict[str, list]:
        """A CSV table of the output (`residuals`, `forces`, `probes`,
        `balance`, `boundaries`, `line_<name>`, `wall_<boundary>`, or a file
        name), as columns; numbers become floats, comment lines (`#`) are
        skipped."""
        p = self.out / (name if name.endswith(".csv") else f"{name}.csv")
        return read_csv(p)

    def line(self, name: str) -> dict[str, list]:
        """The values along `[[output.lines]]` `name`."""
        return self.table(f"line_{name}")

    def record(self) -> list[dict]:
        """The rows of the run record (`run.jsonl`)."""
        with open(self.out / "run.jsonl", encoding="utf-8") as f:
            return [json.loads(line) for line in f if line.strip()]

    # -- fields
    def fields(self, file: str = "fields.vtkhdf", block: str = "internal") -> dict:
        """The cell fields of a VTKHDF file of the run as NumPy arrays
        (needs h5py): `block` is `internal` (the cells) or a boundary name
        (its faces). Adds `points`, the node coordinates."""
        return read_vtkhdf(self.out / file, block)


def _number(s: str):
    try:
        return float(s)
    except ValueError:
        return s


def read_csv(path) -> dict[str, list]:
    """A CSV file of the solver as columns (comment lines skipped)."""
    with open(path, encoding="utf-8") as f:
        rows = [r for r in f if not r.startswith("#") and r.strip()]
    reader = csv.DictReader(rows)
    cols: dict[str, list] = {k: [] for k in (reader.fieldnames or [])}
    for r in reader:
        for k, v in r.items():
            cols[k].append(_number(v))
    return cols


def read_vtkhdf(path, block: str = "internal") -> dict:
    """The cell data of one block of a VTKHDF file as NumPy arrays. For
    `internal` the blocks of the file come back in one piece
    (`foehncraft.vtkhdf.internal`): each cell once, in the order of the cell
    identities, and `points` each once, in the order of the point
    identities."""
    import h5py
    import numpy as np

    from .vtkhdf import internal

    with h5py.File(path, "r") as f:
        if block == "internal":
            g = internal(f)
            return {**g["CellData"], "points": g["Points"]}
        g = f["VTKHDF"][block]
        out = {k: np.asarray(v) for k, v in g["CellData"].items()}
        if "Points" in g:
            out["points"] = np.asarray(g["Points"])
    return out
