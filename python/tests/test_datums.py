"""The tables of the spray datums (SP1b, D267) read as columns."""

import json

import pytest

from foehncraft import Result


def write_output(d):
    (d / "datums.csv").write_text(
        "# format foehncraft-datums 1\n# units: ...\n"
        "step,time,x06.mass_flow,x06.d32,out.mass_flow\n"
        "10,5e-2,1.9e-1,8e-5,1.9e-1\n20,1e-1,1.9e-1,8.1e-5,1.9e-1\n"
    )
    (d / "spray_hist_x06.csv").write_text(
        "# format foehncraft-spray-histogram 1\nstep,time,bin,d_lo,d_hi,droplets,mass\n"
        "10,5e-2,0,1e-6,2e-6,3,1e-12\n"
    )
    (d / "series.json").write_text(
        json.dumps(
            {
                "format": "foehncraft.series",
                "format_version": 1,
                "case": "c",
                "entries": [{"step": 10, "time": 0.05, "files": {"parcels": "spray/p.vtkhdf"}}],
            }
        )
    )
    (d / "result.json").write_text(json.dumps({"status": "finished", "solves": []}))


def test_datum_columns_by_name(tmp_path):
    write_output(tmp_path)
    r = Result(tmp_path)
    x = r.datums("x06")
    assert x["mass_flow"] == [0.19, 0.19]
    assert x["d32"][1] == 8.1e-5
    assert x["step"] == [10, 20]
    assert "out.mass_flow" in r.datums()
    with pytest.raises(KeyError, match="it has out, x06"):
        r.datums("nope")
    assert r.spray_histogram("x06")["droplets"] == [3]
    with pytest.raises(KeyError, match="there are x06"):
        r.spray_histogram("y")
    assert r.series()[0]["files"]["parcels"] == "spray/p.vtkhdf"
