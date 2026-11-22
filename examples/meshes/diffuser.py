"""Structured quadrilateral mesh of the asymmetric plane diffuser of
Buice & Eaton (D45, stage 8) as an SU2 file.

    python3 examples/meshes/diffuser.py OUT.su2 [--level L]

Geometry in units of the inlet channel height H = 1 (Buice & Eaton 1995:
10 degree opening of one wall, area ratio 4.7, so the diffuser is
(4.7 - 1) / tan 10 deg = 20.99 H long). The upper wall is straight at
y = 1, the lower wall is flat at y = 0 over the inlet duct, falls with
10 degrees over the diffuser and is flat at y = -3.7 in the outlet duct.
The lengths of the ducts are our choice (the 1995 brief does not give
them legibly): an inlet duct of 20 H and an outlet duct of 50 H. The
inlet is at x = -20, the diffuser starts at x = 0.

Cells: level 0 has 40 + 42 + 50 intervals in x (inlet duct, diffuser,
outlet duct, the last stretched towards the outlet) and 40 across, graded
with a tanh law towards both walls; every level halves all intervals (the
nodes are nested). Markers: inlet, outlet, upper, lower. The case file
takes `thickness` for the one layer in z.
"""

import math
import sys

H_RATIO = 4.7
ANGLE = math.radians(10.0)
L_IN, L_OUT = 20.0, 50.0
L_DIFF = (H_RATIO - 1.0) / math.tan(ANGLE)
BETA = 2.6          # tanh clustering across the channel
KAPPA = 1.8         # stretching of the outlet duct


def x_nodes(level: int) -> list[float]:
    f = 1 << level
    n_in, n_d, n_o = 40 * f, 42 * f, 50 * f
    xs = [-L_IN + L_IN * i / n_in for i in range(n_in)]
    xs += [L_DIFF * i / n_d for i in range(n_d)]
    xs += [L_DIFF + L_OUT * (math.exp(KAPPA * i / n_o) - 1.0) / (math.exp(KAPPA) - 1.0) for i in range(n_o + 1)]
    return xs


def y_low(x: float) -> float:
    return -min(max(x, 0.0), L_DIFF) * math.tan(ANGLE)


def g(s: float) -> float:
    return 0.5 * (1.0 + math.tanh(BETA * (2.0 * s - 1.0)) / math.tanh(BETA))


def write(path: str, level: int) -> None:
    xs = x_nodes(level)
    nj = 40 << level
    ni = len(xs) - 1
    pts = [(x, y_low(x) + (1.0 - y_low(x)) * g(j / nj)) for x in xs for j in range(nj + 1)]
    idx = lambda i, j: i * (nj + 1) + j  # noqa: E731
    out = ["NDIME= 2", f"NELEM= {ni * nj}"]
    e = 0
    for i in range(ni):
        for j in range(nj):
            out.append(f"9 {idx(i, j)} {idx(i + 1, j)} {idx(i + 1, j + 1)} {idx(i, j + 1)} {e}")
            e += 1
    out.append(f"NPOIN= {len(pts)}")
    out += [f"{x!r} {y!r} {k}" for k, (x, y) in enumerate(pts)]
    out.append("NMARK= 4")
    for name, lines in [
        ("inlet", [(idx(0, j), idx(0, j + 1)) for j in range(nj)]),
        ("outlet", [(idx(ni, j + 1), idx(ni, j)) for j in range(nj)]),
        ("upper", [(idx(i, nj), idx(i + 1, nj)) for i in range(ni)]),
        ("lower", [(idx(i + 1, 0), idx(i, 0)) for i in range(ni)]),
    ]:
        out.append(f"MARKER_TAG= {name}")
        out.append(f"MARKER_ELEMS= {len(lines)}")
        out += [f"3 {a} {b}" for a, b in lines]
    with open(path, "w") as f:
        f.write("\n".join(out) + "\n")


if __name__ == "__main__":
    lvl = int(sys.argv[sys.argv.index("--level") + 1]) if "--level" in sys.argv else 0
    write(sys.argv[1], lvl)
