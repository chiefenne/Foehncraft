"""Writes examples/cases/structure/reed_stop.msh: a curved reed stop (the
guard above a reed that limits its lift) as twenty-node hexahedra in the
Gmsh .msh 4.1 ASCII format, second order, with physical groups.

The stop is a strip bent on a circle: 25 mm long along its centre line,
radius 40 mm, 8 mm wide, 1.5 mm thick, rising away from the reed. The
mesh is mapped from a block, so the edge nodes lie on the arcs. Physical
groups: the volume "stop" and the surfaces "mount" (the clamped end),
"underside" (where the reed lands), "top" and "tip".

No Gmsh installation needed. Run from the repository root:

    python3 -m tools.mesh.reed_stop_msh
"""
import math
import os

OUT = "examples/cases/structure/reed_stop.msh"
LENGTH, RADIUS, WIDTH, THICK = 25e-3, 40e-3, 8e-3, 1.5e-3
NS, NW, NT = 20, 4, 2  # cells along, across, through the thickness

# the twenty-node hexahedron: corners (VTK order), then the edges
CORNERS = [(-1, -1, -1), (1, -1, -1), (1, 1, -1), (-1, 1, -1),
           (-1, -1, 1), (1, -1, 1), (1, 1, 1), (-1, 1, 1)]
EDGES = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4),
         (0, 4), (1, 5), (2, 6), (3, 7)]
# Gmsh node of each of our nodes (crates/foehncraft-structure/src/solid.rs)
GMSH_HEX20 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 13, 9, 16, 18, 19, 17, 10, 12, 14, 15]
# faces: corners then edge nodes (any orientation: the reader orients them)
FACES = {
    "kmin": [0, 3, 2, 1, 11, 10, 9, 8],
    "kmax": [4, 5, 6, 7, 12, 13, 14, 15],
    "imin": [0, 4, 7, 3, 16, 15, 19, 11],
    "imax": [1, 2, 6, 5, 9, 18, 13, 17],
}
NATURAL = CORNERS + [tuple((CORNERS[a][d] + CORNERS[b][d]) / 2 for d in range(3)) for a, b in EDGES]


def position(i, j, k):
    """The point of lattice node (i, j, k): i along the arc, j across, k
    from the underside up."""
    theta = (LENGTH / RADIUS) * i / (2 * NS)
    y = WIDTH * (j / (2 * NW) - 0.5)
    r = RADIUS + THICK / 2 - THICK * k / (2 * NT)
    return (r * math.sin(theta), y, RADIUS - r * math.cos(theta))


def main():
    ids = {}
    nodes = []

    def node(p):
        if p not in ids:
            ids[p] = len(nodes) + 1
            nodes.append(position(*p))
        return ids[p]

    elements = []
    faces = {"mount": [], "underside": [], "top": [], "tip": []}
    side = {"imin": ("mount", 0), "imax": ("tip", NS - 1), "kmin": ("underside", 0), "kmax": ("top", NT - 1)}
    for k in range(NT):
        for j in range(NW):
            for i in range(NS):
                ours = [node((2 * i + int(n[0] + 1), 2 * j + int(n[1] + 1), 2 * k + int(n[2] + 1)))
                        for n in NATURAL]
                gmsh = [0] * 20
                for a, g in enumerate(GMSH_HEX20):
                    gmsh[g] = ours[a]
                elements.append(gmsh)
                for f, (name, at) in side.items():
                    cell = {"imin": i, "imax": i, "kmin": k, "kmax": k}[f]
                    if cell == at:
                        faces[name].append([ours[a] for a in FACES[f]])
    names = list(faces)
    out = ["$MeshFormat", "4.1 0 8", "$EndMeshFormat", "$PhysicalNames", str(len(names) + 1)]
    out += [f'2 {t + 1} "{n}"' for t, n in enumerate(names)]
    out += ['3 10 "stop"', "$EndPhysicalNames", "$Entities", f"0 0 {len(names)} 1"]
    lo = [min(p[d] for p in nodes) for d in range(3)]
    hi = [max(p[d] for p in nodes) for d in range(3)]
    box = " ".join(f"{v:.9g}" for v in lo + hi)
    out += [f"{t + 1} {box} 1 {t + 1} 0" for t in range(len(names))]
    out += [f"1 {box} 1 10 0", "$EndEntities"]
    n = len(nodes)
    out += ["$Nodes", f"1 {n} 1 {n}", f"3 1 0 {n}"]
    out += [str(t) for t in range(1, n + 1)]
    out += [f"{p[0]:.12g} {p[1]:.12g} {p[2]:.12g}" for p in nodes]
    out += ["$EndNodes"]
    total = len(elements) + sum(len(v) for v in faces.values())
    out += ["$Elements", f"{len(names) + 1} {total} 1 {total}"]
    tag = 1
    for t, name in enumerate(names):
        out.append(f"2 {t + 1} 16 {len(faces[name])}")
        for f in faces[name]:
            out.append(f"{tag} " + " ".join(map(str, f)))
            tag += 1
    out.append(f"3 1 17 {len(elements)}")
    for e in elements:
        out.append(f"{tag} " + " ".join(map(str, e)))
        tag += 1
    out.append("$EndElements")
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", encoding="utf-8") as f:
        f.write("\n".join(out) + "\n")
    print(f"wrote {OUT}: {n} nodes, {len(elements)} twenty-node hexahedra")


if __name__ == "__main__":
    main()
