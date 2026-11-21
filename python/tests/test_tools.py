"""The tool registry from Python (D196): every tool once through `call`,
and the catalogue."""

import shutil

from conftest import CASES, needs_native

import foehncraft as fc

ROOT = CASES.parents[2]


@needs_native
def test_the_catalogue_lists_every_tool():
    cat = fc.commands()
    assert cat["format"] == "foehn.tool-catalogue"
    ids = [t["id"] for t in cat["tools"]]
    assert len(ids) == len(set(ids))
    for t in cat["tools"]:
        assert t["input_schema"]["type"] == "object", t["id"]
        assert t["action"] in ("read", "write-files", "start-run"), t["id"]


@needs_native
def test_every_tool_through_call(tmp_path):
    case = tmp_path / "couette.toml"
    shutil.copy(CASES / "couette_hex.toml", case)
    out = tmp_path / "run"
    inputs = {
        "system.info": {},
        "system.backends": {},
        "system.commands": {},
        "case.schema": {},
        "case.check": {"case": case},
        "case.resolve": {"case": case, "set": ["fluid.viscosity=0.02"]},
        "case.estimate": {"case": case, "out": out},
        "case.upgrade": {"files": [case]},
        "diag.codes": {},
        "materials.list": {"temperature": 320.0},
        "materials.show": {"name": "water"},
        "expression.eval": {"source": "2*x + t", "points": [[1.0, 0.0, 0.0, 0.5]]},
        "mesh.process": {
            "input": ROOT / "examples/cases/mesh/hex_2d.toml",
            "out": tmp_path / "mesh",
        },
        "flow.run": {"case": case, "out": out},
        "result.read": {"dir": out},
        "post.derive": {"case": case, "out": out, "fields": ["vorticity"]},
        "study.scalar": {
            "case": ROOT / "examples/cases/scalar/mms_poly_2d.toml",
            "out": tmp_path / "scalar",
        },
        "study.poisson": {
            "case": ROOT / "examples/cases/poisson/poisson_hex_box.toml",
            "levels": [0],
            "solvers": ["amg"],
            "out": tmp_path / "poisson",
        },
        "study.g1": {
            "case": ROOT / "examples/cases/bench/g1_ci_small.toml",
            "out": tmp_path / "g1",
        },
    }
    listed = [t["id"] for t in fc.commands()["tools"]]
    assert sorted(listed) == sorted(inputs), "a tool of the catalogue has no call here"
    results = {}
    for tool in listed:
        r = fc.call(tool, **inputs[tool])
        assert r["format"] == "foehn.tool-result"
        assert r["status"] == "finished", (tool, r["diagnostics"])
        results[tool] = r["result"]
    assert results["expression.eval"]["values"] == [2.5]
    assert "viscosity = 0.02" in results["case.resolve"]["toml"]
    assert results["case.check"]["ok"]
    assert results["flow.run"]["format"] == "foehncraft.result"
    for k in ("format", "status", "exit_code", "case", "files"):
        assert results["result.read"][k] == results["flow.run"][k], k
    assert results["case.estimate"]["cells"] == results["case.check"]["estimate"]["cells"]


@needs_native
def test_wrappers_match_call(tmp_path):
    from foehncraft import _native

    case = CASES / "couette_hex.toml"
    import json

    assert json.loads(_native.check(str(case))) == fc.call("case.check", case=case)["result"]
    assert (
        json.loads(_native.material_database("air", 300.0))
        == (fc.call("materials.show", name="air", temperature=300.0)["result"])
    )


@needs_native
def test_a_bad_call_is_a_diagnostic():
    r = fc.call("flow.rnu")
    assert r["status"] == "error"
    assert r["diagnostics"][0]["code"] == "TOOL001"
    r = fc.call("case.check", cas="x.toml")
    assert r["diagnostics"][0]["code"] == "TOOL002"
