#!/usr/bin/env python3
"""Developed inlet profile of the diffuser case (D45, 7a).

Runs the periodic channel `diffuser_channel_profile.toml` (SST, Re_H =
U_c H / nu = 20 000 on the centreline velocity), reads the profile from
the centre plane to the wall and writes the table `diffuser_inlet_profile.csv`
(y, u_n, k, omega over the full height 0 to H, mirrored about the centre)
that the inlet of `diffuser_sst.toml` reads.

    python3 examples/studies/diffuser_inlet_profile.py [--set mass_flow=0.0893]

The mass flow of the periodic pair sets the bulk velocity; 0.0893 m^2/s
(U_b = 0.893 m/s) gives U_c = 1.000 m/s. Change it, and run this script
again, if the viscosity or the model changes.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))

from foehncraft.profiles import inlet_table, write  # noqa: E402
from foehncraft.results import read_csv  # noqa: E402

CASE = ROOT / "examples/cases/flow/diffuser_channel_profile.toml"
OUT = ROOT / "examples/cases/flow/diffuser_inlet_profile.csv"
H = 1.0


def main() -> int:
    run = ROOT / "output/diffuser_channel_profile"
    cmd = [str(ROOT / "target/release/foehn"), "flow", "--out", str(run), str(CASE)] + sys.argv[1:]
    subprocess.run(cmd, check=False, stdout=subprocess.DEVNULL)
    line = run / "line_profile.csv"
    cols = read_csv(line)
    u_c = float(cols["u"][0])
    half = inlet_table(cols)
    r, u, k, om = half["r"], half["u_n"], half["k"], half["omega"]
    # mirror about the centre plane y = H/2: wall row first
    # (a sample on the centre plane itself is taken once)
    n0 = 1 if r[0] == 0.0 else 0
    ys = [H / 2 - v for v in reversed(r)] + [H / 2 + v for v in r[n0:]]
    table = {
        "y": ys,
        "u_n": list(reversed(u)) + u[n0:],
        "k": list(reversed(k)) + k[n0:],
        "omega": list(reversed(om)) + om[n0:],
    }
    write(table, OUT, header=f"developed channel profile, SST, U_c = {u_c:.4f} m/s, H = {H} m (diffuser_inlet_profile.py)")
    print(f"{OUT.name}: {len(ys)} rows, centreline U_c = {u_c:.4f} m/s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
