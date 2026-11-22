"""Structured quadrilateral mesh of the two-dimensional U-duct of Monson,
Seegmiller & McConnaughey (1990, AIAA 90-1484; D45, stage 8) as an SU2 file.

    python3 examples/meshes/ubend.py OUT.su2 [--level L] [--sharp]

Geometry in units of the duct height H = 1: a 180 degree turn of the
centreline radius Rc = H (Rc/H = 1, inner wall radius 0.5 H, outer wall
radius 1.5 H) about the origin. The inlet leg lies at -1.5 < y < -0.5 and
carries the flow in +x from x = -L_IN to the bend at x = 0; the outlet leg
lies at 0.5 < y < 1.5 and carries it back in -x from the bend to the outlet
at x = -L_OUT. The lengths of the legs are our choice: 6 H upstream
(the inlet takes a developed channel profile) and 20 H downstream (the
experiment has about 14 H; the outlet pressure is kept further away from
the separation behind the bend).

Rounded (default, the experiment): one block that follows the centreline,
i across (j = 0 on the inner wall, j = NJ on the outer one, tanh clustering
towards both walls) and i along: N_IN cells in the inlet leg, stretched
towards the bend, N_BEND cells uniform in angle over the bend, N_OUT cells
in the outlet leg, stretched from the bend towards the outlet. The cells of
the legs next to the bend have the centreline length of the bend cells.

--sharp: the same duct with sharp edges, a divider of thickness H with a
square end at x = 0 and a square end wall at x = H (a turn gap of H);
Cartesian cells on a tensor grid without the cells of the divider, so the
two convex corners of the divider end and the two concave corners of the
end wall are kinks. The cells of the turn gap are clustered towards the
divider end and the end wall like the cells across the duct, and the legs
start from that spacing, so the cells at the kinks are square
(docs/pitfalls.md, diffuser). The inlet is at the grid line next to
x = -L_IN.

Every level halves all intervals (the nodes are nested). Markers: inlet,
outlet, inner, outer (sharp: the divider is `inner`, the outer walls and
the end wall are `outer`). The case file takes `thickness` for the one
layer in z.
"""

import math
import sys

L_IN, L_OUT = 6.0, 20.0
RI, RO = 0.5, 1.5
FIRST = 8e-4        # wall-cell height on level 0 (y+ about 4 at Re_b 1e5; 1 on level 2)
N_IN, N_BEND, N_OUT, N_J = 24, 96, 64, 48
# --sharp: turn gap of H with N_GAP cells clustered towards the divider end
# and the end wall like the cells across the duct
GAP, N_GAP = 1.0, 48
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


def wall_beta(n: int, first: float) -> float:
    """The tanh clustering of g() whose first interval over a unit height is `first`."""
    lo, hi = 0.1, 10.0
    for _ in range(200):
        b = 0.5 * (lo + hi)
        h = 0.5 * (1.0 + math.tanh(b * (2.0 / n - 1.0)) / math.tanh(b))
        lo, hi = (b, hi) if h > first else (lo, b)
    return 0.5 * (lo + hi)


def across(n: int) -> list[float]:
    """n intervals over [0, 1], tanh clustering towards both ends; the
    clustering is fixed on level 0 and the levels are nested."""
    b = wall_beta(N_J, FIRST)
    return [0.5 * (1.0 + math.tanh(b * (2.0 * j / n - 1.0)) / math.tanh(b)) for j in range(n + 1)]


def leg(n: int, length: float, first: float) -> list[float]:
    """Distances 0 .. length from the bend, the first interval `first`."""
    k = stretch(first, n, length)
    return [length * math.expm1(k * i / n) / math.expm1(k) for i in range(n + 1)]


def rounded(level: int):
    f = 1 << level
    n_in, n_b, n_o, nj = N_IN * f, N_BEND * f, N_OUT * f, N_J * f
    ds = math.pi / n_b                       # centreline length of a bend cell (Rc = 1)
    eta = across(nj)
    s_in = [-d for d in reversed(leg(n_in, L_IN, ds))]     # -L_IN .. 0
    phi = [-0.5 * math.pi + math.pi * i / n_b for i in range(n_b + 1)]
    s_out = leg(n_o, L_OUT, ds)                            # 0 .. L_OUT
    cols = []
    for x in s_in[:-1]:
        cols.append([(x, -(RI + e)) for e in eta])
    for p in phi:
        cols.append([((RI + e) * math.cos(p), (RI + e) * math.sin(p)) for e in eta])
    for d in s_out[1:]:
        cols.append([(-d, RI + e) for e in eta])
    ni = len(cols) - 1
    pts = [pt for c in cols for pt in c]
    idx = lambda i, j: i * (nj + 1) + j  # noqa: E731
    quads = []
    for i in range(ni):
        for j in range(nj):
            # inlet leg: i runs in +x, j in -y; keep counter-clockwise
            quads.append((idx(i, j), idx(i, j + 1), idx(i + 1, j + 1), idx(i + 1, j)))
    marks = {
        "inlet": [(idx(0, j + 1), idx(0, j)) for j in range(nj)],
        "outlet": [(idx(ni, j), idx(ni, j + 1)) for j in range(nj)],
        "inner": [(idx(i + 1, 0), idx(i, 0)) for i in range(ni)],
        "outer": [(idx(i, nj), idx(i + 1, nj)) for i in range(ni)],
    }
    return pts, quads, marks


def sharp(level: int):
    f = 1 << level
    n_o, n_g, nj = N_OUT * f, N_GAP * f, N_J * f
    xg = [GAP * e for e in across(n_g)]
    dx0 = xg[1]
    # the legs: one x grid for both (the outlet leg is the longer one)
    xs_leg = [-d for d in reversed(leg(n_o, L_OUT, dx0))]   # -L_OUT .. 0
    xs = xs_leg[:-1] + xg
    eta = across(nj)
    ys = [-RO + e for e in eta] + [-RI + e for e in eta[1:]] + [RI + e for e in eta[1:]]
    nx, ny = len(xs), len(ys)
    jlo, jhi = nj, 2 * nj                     # rows of y = -0.5 and y = 0.5
    i0 = len(xs_leg) - 1                      # column of x = 0
    iin = min(range(nx), key=lambda i: abs(xs[i] + L_IN))   # inlet column near x = -L_IN
    def cell(i, j):
        if i < iin and j < jlo:
            return False                      # inlet leg shorter than the outlet leg
        return not (i < i0 and jlo <= j < jhi)    # the divider
    used = {}
    pts = []
    def node(i, j):
        k = i * ny + j
        if k not in used:
            used[k] = len(pts)
            pts.append((xs[i], ys[j]))
        return used[k]
    quads, cid = [], {}
    for i in range(nx - 1):
        for j in range(ny - 1):
            if cell(i, j):
                cid[(i, j)] = len(quads)
                quads.append((node(i, j), node(i + 1, j), node(i + 1, j + 1), node(i, j + 1)))
    marks = {"inlet": [], "outlet": [], "inner": [], "outer": []}
    for (i, j) in cid:
        if j == 0:                                    # y = -1.5
            marks["outer"].append((node(i + 1, j), node(i, j)))
        if j == ny - 2:                               # y = 1.5
            marks["outer"].append((node(i, j + 1), node(i + 1, j + 1)))
        if i == nx - 2:                               # end wall x = GAP
            marks["outer"].append((node(i + 1, j), node(i + 1, j + 1)))
        if i == 0 or not cell(i - 1, j):              # inlet, outlet or divider end
            name = "inlet" if j < jlo else ("outlet" if j >= jhi else "inner")
            marks[name].append((node(i, j + 1), node(i, j)))
        if j > 0 and not cell(i, j - 1):              # divider, upper side
            marks["inner"].append((node(i + 1, j), node(i, j)))
        if j < ny - 2 and not cell(i, j + 1):         # divider, lower side
            marks["inner"].append((node(i, j + 1), node(i + 1, j + 1)))
    return pts, quads, marks


def write(path: str, level: int) -> None:
    pts, quads, marks = sharp(level) if SHARP else rounded(level)
    out = ["NDIME= 2", f"NELEM= {len(quads)}"]
    out += [f"9 {a} {b} {c} {d} {e}" for e, (a, b, c, d) in enumerate(quads)]
    out.append(f"NPOIN= {len(pts)}")
    out += [f"{x!r} {y!r} {k}" for k, (x, y) in enumerate(pts)]
    out.append(f"NMARK= {len(marks)}")
    for name, lines in marks.items():
        out.append(f"MARKER_TAG= {name}")
        out.append(f"MARKER_ELEMS= {len(lines)}")
        out += [f"3 {a} {b}" for a, b in lines]
    with open(path, "w") as fh:
        fh.write("\n".join(out) + "\n")


if __name__ == "__main__":
    lvl = int(sys.argv[sys.argv.index("--level") + 1]) if "--level" in sys.argv else 0
    SHARP = "--sharp" in sys.argv
    write(sys.argv[1], lvl)
