#!/usr/bin/env python3
"""A thin sector of a sphere as a chain of hexahedra (SU2 format), for
spherically symmetric problems (the Scriven bubble of stage MP2).

The four side planes go through the centre, so a radial flow is tangent to
them and `symmetry` on the marker `side` is exact; the cells are planar
frustums of a square pyramid of half angle `--angle` about +z. The marker
`centre` closes the chain at the small radius `--r0`, `outer` at `--r1`.

    python3 examples/meshes/sphere_sector.py OUT.su2 --r0 2e-6 --r1 6e-4 --cells 600
"""

import argparse
import math


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("out")
    ap.add_argument("--r0", type=float, default=2.0e-6, help="inner radius (m)")
    ap.add_argument("--r1", type=float, default=6.0e-4, help="outer radius (m)")
    ap.add_argument("--cells", type=int, default=600, help="radial cells")
    ap.add_argument("--angle", type=float, default=2.0, help="half angle (degrees)")
    a = ap.parse_args()
    t = math.tan(math.radians(a.angle))
    n = math.sqrt(1.0 + 2.0 * t * t)
    # corners counter-clockwise seen from the centre: the normal of the
    # inner quad points outwards (VTK hexahedron order)
    dirs = [(-t / n, -t / n, 1 / n), (t / n, -t / n, 1 / n), (t / n, t / n, 1 / n), (-t / n, t / n, 1 / n)]
    rs = [a.r0 + (a.r1 - a.r0) * i / a.cells for i in range(a.cells + 1)]
    pts = [(r * d[0], r * d[1], r * d[2]) for r in rs for d in dirs]
    hexes = [[4 * i + k for k in range(4)] + [4 * (i + 1) + k for k in range(4)] for i in range(a.cells)]
    sides = []
    for i in range(a.cells):
        for k in range(4):
            j = (k + 1) % 4
            sides.append((4 * i + k, 4 * (i + 1) + k, 4 * (i + 1) + j, 4 * i + j))
    m = 4 * a.cells
    with open(a.out, "w") as f:
        f.write("NDIME= 3\n")
        f.write(f"NELEM= {len(hexes)}\n")
        for i, h in enumerate(hexes):
            f.write("12 " + " ".join(map(str, h)) + f" {i}\n")
        f.write(f"NPOIN= {len(pts)}\n")
        for i, p in enumerate(pts):
            f.write(f"{p[0]:.17e} {p[1]:.17e} {p[2]:.17e} {i}\n")
        f.write("NMARK= 3\n")
        f.write("MARKER_TAG= centre\nMARKER_ELEMS= 1\n9 3 2 1 0\n")
        f.write(f"MARKER_TAG= outer\nMARKER_ELEMS= 1\n9 {m} {m + 1} {m + 2} {m + 3}\n")
        f.write(f"MARKER_TAG= side\nMARKER_ELEMS= {len(sides)}\n")
        for s in sides:
            f.write("9 " + " ".join(map(str, s)) + "\n")


if __name__ == "__main__":
    main()
