"""The STL surface of the blown cylinder (L13 stage 2e, gallery `cut-patches`).

A circular cylinder of radius 0.05 m about (0.2, 0.2), as a closed prism
from z = -0.005 to 0.015 m (it spans the 0.01 m deep 2D mesh), in two
named solids of an ASCII STL file: `slot`, the part of the mantle from 45°
to 75° (the blowing slot), and `wall`, the rest with both caps. 120 facets
around (3°), outward normals.

    python3 -m tools.mesh.blown_cylinder_stl examples/cases/flow/immersed_blown_cylinder.stl
"""
import math
import sys

CX, CY, R = 0.2, 0.2, 0.05
Z0, Z1 = -0.005, 0.015
N = 120
SLOT = (45.0, 75.0)


def facets():
    """(solid, triangle) pairs, outward."""
    p = [(CX + R * math.cos(2 * math.pi * i / N), CY + R * math.sin(2 * math.pi * i / N))
         for i in range(N)]
    out = []
    for i in range(N):
        a, b = p[i], p[(i + 1) % N]
        mid = (i + 0.5) * 360.0 / N
        solid = "slot" if SLOT[0] < mid < SLOT[1] else "wall"
        a0, b0, a1, b1 = (*a, Z0), (*b, Z0), (*a, Z1), (*b, Z1)
        out.append((solid, (a0, b0, b1)))
        out.append((solid, (a0, b1, a1)))
        out.append(("wall", ((CX, CY, Z0), b0, a0)))
        out.append(("wall", ((CX, CY, Z1), a1, b1)))
    return out


def write(path):
    f = facets()
    with open(path, "w") as o:
        for solid in ("wall", "slot"):
            o.write(f"solid {solid}\n")
            for s, t in f:
                if s != solid:
                    continue
                u = [t[1][k] - t[0][k] for k in range(3)]
                v = [t[2][k] - t[0][k] for k in range(3)]
                n = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]
                ln = math.sqrt(sum(x * x for x in n))
                o.write("  facet normal {:.9e} {:.9e} {:.9e}\n    outer loop\n".format(*(x / ln for x in n)))
                for q in t:
                    o.write("      vertex {:.9e} {:.9e} {:.9e}\n".format(*q))
                o.write("    endloop\n  endfacet\n")
            o.write(f"endsolid {solid}\n")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    write(sys.argv[1])
