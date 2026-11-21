"""Runs from Python (D61): in process and through the command line, the
written case file is the run, callbacks see the run record and can stop."""

import json

import h5py
import numpy as np
import pytest
from conftest import CASES, foehn_path, needs_foehn, needs_native

from foehncraft import Case, FoehncraftError, Result, run
from foehncraft.expr import constant, y


def datasets(path):
    """Every dataset of an HDF5 file as bytes."""
    out = {}
    with h5py.File(path, "r") as f:
        f.visititems(
            lambda n, o: (
                out.__setitem__(n, np.asarray(o).tobytes()) if isinstance(o, h5py.Dataset) else None
            )
        )
    return out


def same_fields(a, b):
    da, db = datasets(a / "fields.vtkhdf"), datasets(b / "fields.vtkhdf")
    assert da.keys() == db.keys()
    for k in da:
        assert da[k] == db[k], k


@needs_native
def test_a_case_built_in_python_is_its_case_file(tmp_path):
    # the expression inlet of poiseuille_expression.toml, built in Python
    c = Case.load(CASES / "poiseuille_hex.toml")
    c["verification"] = {"max_error": 2e-3}
    c["constants"] = {"u_mean": 1.0, "h": 1.0}
    um, h = constant("u_mean"), constant("h")
    c["boundary.xmin.velocity"] = [6 * um * (y / h) * (1 - y / h), 0.0, 0.0]
    del c["boundary.xmin.profile"]
    c["output.vtkhdf"] = True
    a = run(c, out=tmp_path / "py")
    assert a.ok, a.diagnostics
    # one mesh: a limitation (D199), the run still finished with exit code 0
    assert a.status == "finished"
    assert [lim["code"] for lim in a.limitations] == ["LIM004"]
    # the same case written by hand as TOML strings
    text = (CASES / "poiseuille_hex.toml").read_text()
    text = text.replace(
        'profile = "analytic"     # the exact parabola',
        'velocity = ["6 * u_mean * (y / h) * (1 - y / h)", 0.0, 0.0]',
    ).replace("levels = [0, 1, 2]\nexpected_order = [1.8, 2.2]\n", "max_error = 2e-3\n")
    text = text.replace("vtkhdf = false", "vtkhdf = true")
    text = "[constants]\nu_mean = 1.0\nh = 1.0\n\n" + text
    hand = tmp_path / "hand.toml"
    hand.write_text(text)
    b = run(hand, out=tmp_path / "hand")
    assert b.ok
    same_fields(a.out, b.out)
    assert a.value("errors.l2") == b.value("errors.l2")


@needs_native
def test_callbacks_see_every_iteration_and_stop_the_run(tmp_path):
    rows = []

    def on_iteration(row):
        rows.append(row)
        return row["iteration"] >= 5

    res = run(CASES / "cavity_re100.toml", out=tmp_path / "stop", on_iteration=on_iteration)
    assert res.status == "stopped" and res.exit_code == 3
    assert [r["iteration"] for r in rows] == [1, 2, 3, 4, 5]
    assert {"momentum", "continuity", "cg"} <= rows[0].keys()
    # the record has the same rows
    it = [r for r in res.record() if r["type"] == "iteration"]
    assert [r["iteration"] for r in it] == [1, 2, 3, 4, 5]


@needs_native
def test_steps_of_a_transient_run(tmp_path):
    steps = []
    res = run(
        CASES / "taylor_green_bdf2.toml",
        out=tmp_path / "tg",
        sets={"time.end_time": 0.02, "verification.time_steps": [], "verification.max_error": 1.0},
        on_step=lambda r: steps.append((r["step"], r["time"])),
    )
    assert res.exit_code in (0, 2), res.diagnostics
    assert steps and steps[-1][1] == pytest.approx(0.02)


@needs_native
def test_an_exception_in_a_callback_is_raised(tmp_path):
    def bad(row):
        raise ValueError("from the callback")

    with pytest.raises(ValueError, match="from the callback"):
        run(CASES / "couette_hex.toml", out=tmp_path / "bad", on_iteration=bad)


@needs_native
def test_errors_of_the_case_carry_codes(tmp_path):
    c = Case.load(CASES / "couette_hex.toml")
    c["fluid.material"] = "watr"
    with pytest.raises(FoehncraftError) as ex:
        run(c, out=tmp_path / "e")
    assert ex.value.code == "MAT001"


@needs_native
@needs_foehn
def test_the_command_line_gives_the_same_run(tmp_path):
    sets = ["output.vtkhdf=true"]
    a = run(CASES / "couette_poly.toml", out=tmp_path / "native", sets=sets, backend="native")
    rows = []
    b = run(
        CASES / "couette_poly.toml",
        out=tmp_path / "cli",
        sets=sets,
        backend="cli",
        foehn=foehn_path(),
        on_iteration=rows.append,
    )
    assert a.ok and b.ok
    same_fields(a.out, b.out)
    assert a.summary["solves"] == b.summary["solves"]
    assert rows, "the command line run record was followed"
    assert json.loads((b.out / "result.json").read_text())["status"] == "finished"


@needs_native
def test_result_tables(tmp_path):
    res = run(CASES / "couette_hex.toml", out=tmp_path / "t")
    r = Result.load(res.out)
    hist = r.table("residuals")
    assert hist["iteration"][0] == 1.0 and len(hist["continuity"]) == len(hist["iteration"])
    assert "errors.l2" in r.flat()
    assert r.files["result.json"] if "result.json" in r.files else True


@needs_native
def test_uncertainty_of_a_sweep(tmp_path):
    # three meshes of the Couette flow with a viscosity gradient (D67)
    res = run(CASES / "couette_viscosity.toml", out=tmp_path / "gci")
    assert res.ok
    e = res.uncertainty["boundaries.ymin.force_x"]
    assert e["verdict"] == "monotonic"
    assert 1.5 < e["order"] < 2.5
    exact = 0.1 * 4.0 / np.log1p(4.0) * 0.2
    assert abs(e["values"][0] - exact) <= e["band"]
    assert abs(e["extrapolated"] - exact) < 0.1 * abs(e["values"][0] - exact)
