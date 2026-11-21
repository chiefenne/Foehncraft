"""The tool registry from Python (D196): every engine function by its id.

    import foehncraft as fc

    cat = fc.commands()                      # the catalogue: ids, schemas, action kinds
    r = fc.call("case.check", case="examples/cases/flow/cavity_re100.toml",
                set=["fluid.viscosity=1e-3"])
    r["status"], r["result"]["ok"]

`call` returns the result envelope (`foehn.tool-result`) as a dict: the
tool, its version, `status`, `exit_code`, `result` (the tool's output) and
`diagnostics` (the error of a call that failed; it does not raise). The
input fields are those of the tool's input schema; paths may be `Path`
objects. The same tools are what `foehn call` and the command line reach.
With the native module the call runs in process; without it, through
`foehn call`.
"""

from __future__ import annotations

import json
import subprocess
from collections.abc import Callable
from typing import Any

from .run import _backend, find_foehn


def _json(fields: dict[str, Any]) -> str:
    return json.dumps(fields, default=str)


def call(
    tool: str,
    /,
    *,
    observer: Callable[[dict], Any] | None = None,
    log: Callable[[str], Any] | None = None,
    backend: str = "auto",
    foehn: str | None = None,
    **fields: Any,
) -> dict:
    """Runs the tool `tool` with the input `fields` and returns the envelope.

    `observer(row)` sees every row of the run record of a run (a true value
    stops it), `log(line)` the human text; both only with the native module.
    """
    if _backend(backend) == "native":
        from . import _native

        obs = None
        if observer is not None:

            def obs(text: str) -> bool:
                return bool(observer(json.loads(text)))

        return json.loads(_native.call(tool, _json(fields), obs, log))
    p = subprocess.run(
        [find_foehn(foehn), "--quiet", "call", tool, "--input", "-"],
        input=_json(fields),
        capture_output=True,
        text=True,
        check=False,
    )
    return json.loads(p.stdout)


def commands(*, backend: str = "auto", foehn: str | None = None) -> dict:
    """The catalogue of every tool (`foehn.tool-catalogue`)."""
    if _backend(backend) == "native":
        from . import _native

        return json.loads(_native.commands())
    p = subprocess.run([find_foehn(foehn), "commands"], capture_output=True, text=True, check=True)
    return json.loads(p.stdout)
