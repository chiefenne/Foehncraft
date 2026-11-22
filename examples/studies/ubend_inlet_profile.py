#!/usr/bin/env python3
"""Developed inlet profile of the U-bend case (D45, stage 8).

Runs the periodic channel `ubend_channel_profile.toml` (SST, Re_b =
U_b H / nu = 1e5), reads the profile from the centre plane to the wall and
writes the table `ubend_inlet_profile.csv` (y, u_n, k, omega over the full
height of the inlet leg, -1.5 H to -0.5 H, mirrored about the centre) that
the inlet of `ubend_sst.toml` and `ubend_sst_sharp.toml` reads.

    python3 examples/studies/ubend_inlet_profile.py
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))

from foehncraft.profiles import inlet_table, write  # noqa: E402
from foehncraft.results import read_csv  # noqa: E402

CASE = ROOT / "examples/cases/flow/ubend_channel_profile.toml"
OUT = ROOT / "examples/cases/flow/ubend_inlet_profile.csv"
H = 1.0
Y0 = -1.5   # the outer wall of the inlet leg


def main() -> int:
    run = ROOT / "output/ubend_channel_profile"
    cmd = [str(ROOT / "target/release/foehn"), "flow", "--out", str(run), str(CASE)] + sys.argv[1:]
    subprocess.run(cmd, check=False, stdout=subprocess.DEVNULL)
    cols = read_csv(run / "line_profile.csv")
    u_c = float(cols["u"][0])
    half = inlet_table(cols)
    r, u, k, om = half["r"], half["u_n"], half["k"], half["omega"]
    n0 = 1 if r[0] == 0.0 else 0
    ys = [Y0 + H / 2 - v for v in reversed(r)] + [Y0 + H / 2 + v for v in r[n0:]]
    table = {
        "y": ys,
        "u_n": list(reversed(u)) + u[n0:],
        "k": list(reversed(k)) + k[n0:],
        "omega": list(reversed(om)) + om[n0:],
    }
    write(table, OUT, header=f"developed channel profile, SST, Re_b 1e5, U_b = 1 m/s, U_c = {u_c:.4f} m/s, H = {H} m (ubend_inlet_profile.py)")
    print(f"{OUT.name}: {len(ys)} rows, centreline U_c = {u_c:.4f} m/s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
