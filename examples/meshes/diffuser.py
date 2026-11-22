"""Structured quadrilateral mesh of the asymmetric plane diffuser of
Buice & Eaton (D45, stage 8) as an SU2 file.

    python3 examples/meshes/diffuser.py OUT.su2 [--level L] [--sharp]

Geometry in units of the inlet channel height H = 1 (Buice & Eaton 1995;
Obi et al. 1993: 10 degree opening of one wall, area ratio 4.7, so the
diffuser is (4.7 - 1) / tan 10 deg = 20.99 H long between the
intersections of the tangents). The upper wall is straight at y = 1, the
lower wall is flat at y = 0 over the inlet duct, falls with 10 degrees
over the diffuser and is flat at y = -3.7 in the outlet duct. Both corners
of the lower wall are rounded with the radius 4.3 H of the experiment
(8.6 half-heights, Kaltenbach 1994, CTR Annual Research Briefs, p. 176):
a sharp corner makes the least-squares velocity gradient of the thin wall
cells wrong by O(1) on every mesh, and the linear-upwind convection
extrapolates that error until the run diverges (docs/pitfalls.md). The
lengths of the ducts are our choice: an inlet duct of 20 H and an outlet
duct of 50 H. The inlet is at x = -20, the tangents meet at x = 0 and
x = 20.99.

Cells: level 0 has 40 + 84 + 50 intervals in x (inlet duct, stretched
towards the diffuser; diffuser, uniform at 0.25 H, so that the wall cells,
trapezoids that grow with the channel, keep a face skewness below 0.5 and
the D50 remediation stays off on every level; outlet duct, stretched from
0.25 H towards the outlet) and 40 across, graded with a tanh law towards
both walls; every level halves all intervals (the nodes are nested).
Markers: inlet, outlet, upper, lower. The case file takes `thickness` for
the one layer in z.
"""

import math
import sys

H_RATIO = 4.7
ANGLE = math.radians(10.0)
TAN = math.tan(ANGLE)
L_IN, L_OUT = 20.0, 50.0
L_DIFF = (H_RATIO - 1.0) / TAN
RADIUS = 4.3        # corner radius (Kaltenbach 1994: 8.6 delta, H = 2 delta)
BETA = 2.6          # tanh clustering across the channel
N_IN, N_DIFF, N_OUT, N_Y = 40, 84, 50, 40
# --sharp: corners without radius; the diffuser intervals shrink towards
# both corners, x = L (s - A sin(2 pi s)/(2 pi)), from 0.33 H in the middle
# to 0.0175 H at the corners on level 0 (wall cells of aspect ratio about 11
# at the kinks instead of 150)
N_DIFF_SHARP, A_SHARP = 120, 0.9
SHARP = False


def stretch(first: float, n: int, length: float) -> float:
    """The k of x = length (e^{k s} - 1)/(e^k - 1), s = i/n, whose first
    interval is `first` (bisection; k < 0 when `first` is above the mean)."""
    def f(k: float) -> float:
        return length * math.expm1(k / n) / math.expm1(k) - first
    lo, hi = -30.0, 30.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if abs(mid) < 1e-12:
            mid = 1e-9
        if f(mid) > 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def x_nodes(level: int) -> list[float]:
    f = 1 << level
    n_in, n_d, n_o = N_IN * f, (N_DIFF_SHARP if SHARP else N_DIFF) * f, N_OUT * f
    dx = L_DIFF / n_d * (1.0 - A_SHARP if SHARP else 1.0)
    k_in = stretch(dx, n_in, L_IN)
    xs = [-L_IN * math.expm1(k_in * (n_in - i) / n_in) / math.expm1(k_in) for i in range(n_in)]
    if SHARP:
        xs += [L_DIFF * (i / n_d - A_SHARP * math.sin(2.0 * math.pi * i / n_d) / (2.0 * math.pi)) for i in range(n_d)]
    else:
        xs += [dx * i for i in range(n_d)]
    k_o = stretch(dx, n_o, L_OUT)
    xs += [L_DIFF + L_OUT * math.expm1(k_o * i / n_o) / math.expm1(k_o) for i in range(n_o + 1)]
    return xs


def y_low(x: float) -> float:
    """The lower wall: straight pieces joined by two arcs of RADIUS, tangent
    to both pieces (the tangents meet at x = 0 and x = L_DIFF)."""
    if SHARP:
        return -min(max(x, 0.0), L_DIFF) * TAN
    t = RADIUS * math.tan(0.5 * ANGLE)
    s = RADIUS * math.sin(ANGLE)
    if -t <= x <= s - t:
        return -RADIUS + math.sqrt(RADIUS * RADIUS - (x + t) ** 2)
    if L_DIFF + t - s <= x <= L_DIFF + t:
        return -3.7 + RADIUS - math.sqrt(RADIUS * RADIUS - (x - L_DIFF - t) ** 2)
    return -min(max(x, 0.0), L_DIFF) * TAN


def g(s: float) -> float:
    return 0.5 * (1.0 + math.tanh(BETA * (2.0 * s - 1.0)) / math.tanh(BETA))


def write(path: str, level: int) -> None:
    xs = x_nodes(level)
    nj = N_Y << level
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
    SHARP = "--sharp" in sys.argv
    write(sys.argv[1], lvl)
