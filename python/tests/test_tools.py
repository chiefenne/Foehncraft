"""The tool registry from Python (D196): every tool once through `call`,
and the catalogue."""

import json
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
        assert t["action"] in ("read", "write-files", "start-run", "control"), t["id"]


@needs_native
def test_every_tool_through_call(tmp_path):
    case = tmp_path / "couette.toml"
    shutil.copy(CASES / "couette_hex.toml", case)
    out = tmp_path / "run"
    # a short spray run with datums for the post.* tools of SP1b
    flow = ROOT / "examples/cases/flow"
    for f in ("spray_datums.toml", "spray_datums_injectors.toml", "spray_datums_section.stl"):
        shutil.copy(flow / f, tmp_path / f)
    spray = tmp_path / "spray_datums.toml"
    spray_out = tmp_path / "spray"
    r = fc.call("flow.run", case=spray, out=spray_out, end_time=0.05)
    assert r["status"] == "finished", r["diagnostics"]
    parcels = (
        spray_out
        / json.loads((spray_out / "series.json").read_text())["entries"][-1]["files"]["parcels"]
    )
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
        "result.page": {"dir": out, "series": "residuals", "limit": 5},
        # live control (L19): no run listens in a finished run's folder
        "run.status": {"dir": out},
        "run.stop": {"dir": out},
        "run.set": {"dir": out, "set": {"time.time_step": 0.01}},
        "post.derive": {"case": case, "out": out, "fields": ["vorticity"]},
        "post.datums": {"case": spray, "out": spray_out, "datums": ["x06", "outlet"]},
        "post.droplets": {"case": spray, "file": parcels, "datum": "nozzle_box"},
        "post.penetration": {"case": spray, "file": parcels},
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
        "structure.run": {
            "case": ROOT / "examples/cases/structure/reed_beam.toml",
            "out": tmp_path / "structure",
        },
        "structure.fitclamp": {
            "case": ROOT / "examples/cases/structure/reed_beam_clamp_springs.toml",
            "frequency": 150.0,
        },
        # the 0D reactor and the equilibrium of stage SF2a
        "chemistry.reactor": {
            "case": ROOT / "examples/cases/chemistry/reactor_methane_air_jl4.toml",
            "csv": tmp_path / "reactor.csv",
        },
        "chemistry.equilibrium": {
            "case": ROOT / "examples/cases/chemistry/reactor_methane_air_gri30.toml",
        },
        # the call of flow.run in its run record, run again
        "journal.replay": {"script": out / "run.jsonl"},
    }
    listed = [t["id"] for t in fc.commands()["tools"]]
    assert sorted(listed) == sorted(inputs), "a tool of the catalogue has no call here"
    results = {}
    for tool in listed:
        r = fc.call(tool, **inputs[tool])
        assert r["format"] == "foehn.tool-result"
        if tool.startswith("run."):
            assert r["status"] == "error", tool
            assert r["diagnostics"][0]["code"] == "TOOL007", tool
            continue
        assert r["status"] == "finished", (tool, r["diagnostics"])
        # the limitations of the output in the envelope (D199)
        assert r["limitations"] == (r["result"] or {}).get("limitations", []), tool
        results[tool] = r["result"]
    assert results["expression.eval"]["values"] == [2.5]
    assert "viscosity = 0.02" in results["case.resolve"]["toml"]
    assert results["case.check"]["ok"]
    assert results["flow.run"]["format"] == "foehncraft.result"
    for k in ("format", "status", "exit_code", "case", "files"):
        assert results["result.read"][k] == results["flow.run"][k], k
    assert results["case.estimate"]["cells"] == results["case.check"]["estimate"]["cells"]
    assert {d["name"] for d in results["post.datums"]["datums"]} == {"x06", "outlet"}
    assert results["post.penetration"]["injectors"][0]["name"] == "nozzle"
    assert results["structure.run"]["modal"]["sturm"]["ok"]
    assert abs(results["structure.fitclamp"]["frequency"] - 150.0) < 1e-3
    assert results["chemistry.reactor"]["failures"] == []
    assert 2690.0 < results["chemistry.equilibrium"]["state"]["temperature"] < 2705.0


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


@needs_native
def test_calls_go_into_the_journal_and_replay(tmp_path, monkeypatch):
    journal = tmp_path / "session.jsonl"
    monkeypatch.setenv("FOEHN_JOURNAL", str(journal))
    fc.call("materials.show", name="water", temperature=330.0)
    fc.call("expression.eval", source="2*x", points=[[1.0, 0.0, 0.0, 0.0]])
    fc.call("case.chek")
    monkeypatch.delenv("FOEHN_JOURNAL")
    rows = [json.loads(line) for line in journal.read_text().splitlines()]
    assert rows[0]["format"] == "foehn.journal"
    calls = rows[1:]
    assert [c["tool"] for c in calls] == ["materials.show", "expression.eval", "case.chek"]
    assert [c["seq"] for c in calls] == [1, 2, 3]
    assert calls[0]["via"] == "python"
    # the typed input: defaults filled in
    assert calls[0]["input"]["pressure"] == 101325.0
    assert calls[2]["diagnostic"]["code"] == "TOOL001"
    r = fc.call("journal.replay", script=journal)
    assert r["exit_code"] == 0, r
    assert (r["result"]["replayed"], r["result"]["skipped"]) == (2, 1)
    # nothing journaled after the variable is gone
    assert len(journal.read_text().splitlines()) == 4
