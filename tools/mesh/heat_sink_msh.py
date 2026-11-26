"""Writes examples/cases/structure/heat_sink.msh: the aluminium heat sink of
examples/cases/flow/cht_heat_sink_turbulent.toml (a base 100 mm long and
5 mm thick with three fins 3 mm thick and 20 mm tall, 40 mm apart, 10 mm
deep) as twenty-node hexahedra in the Gmsh .msh 4.1 ASCII format, second
order, with physical groups: the volume "sink" and the surface "bottom"
(the underside of the base, on the chip and the board).

The structural mesh is its own: finer across the fins and with four
elements through the depth, where the flow mesh has one layer of cells.
No Gmsh installation needed. Run from the repository root:

    python3 -m tools.mesh.heat_sink_msh
"""
import os

OUT = "examples/cases/structure/heat_sink.msh"
# breaks of the lattice and elements per interval
X = [0.080, 0.088, 0.091, 0.128, 0.131, 0.168, 0.171, 0.180]
NX = [6, 3, 20, 3, 20, 3, 6]
Y = [0.0, 0.005, 0.025]
NY = [3, 12]
Z = [0.0, 0.01]
NZ = [4]
FINS = {1, 3, 5}  # x intervals that carry a fin

CORNERS = [(-1, -1, -1), (1, -1, -1), (1, 1, -1), (-1, 1, -1),
           (-1, -1, 1), (1, -1, 1), (1, 1, 1), (-1, 1, 1)]
EDGES = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4),
         (0, 4), (1, 5), (2, 6), (3, 7)]
# Gmsh node of each of our nodes (crates/foehncraft-structure/src/solid.rs)
GMSH_HEX20 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 13, 9, 16, 18, 19, 17, 10, 12, 14, 15]
# the face at y = min of an element: corners then edge nodes
YMIN = [0, 1, 5, 4, 8, 17, 12, 16]
NATURAL = CORNERS + [tuple((CORNERS[a][d] + CORNERS[b][d]) / 2 for d in range(3)) for a, b in EDGES]


def lattice(breaks, counts):
    """Coordinates of the half-element lattice (corner and mid-edge nodes)
    and the interval of each element."""
    pts, owner = [breaks[0]], []
    for s, (a, b) in enumerate(zip(breaks, breaks[1:])):
        n = counts[s]
        for e in range(n):
            pts.append(a + (b - a) * (2 * e + 1) / (2 * n))
            pts.append(a + (b - a) * (2 * e + 2) / (2 * n))
            owner.append(s)
    return pts, owner


def main():
    xs, xo = lattice(X, NX)
    ys, yo = lattice(Y, NY)
    zs, _ = lattice(Z, NZ)
    ids, nodes = {}, []

    def node(p):
        if p not in ids:
            ids[p] = len(nodes) + 1
            nodes.append((xs[p[0]], ys[p[1]], zs[p[2]]))
        return ids[p]

    elements, bottom = [], []
    for k in range(len(zs) // 2):
        for j in range(len(ys) // 2):
            for i in range(len(xs) // 2):
                if yo[j] == 1 and xo[i] not in FINS:
                    continue
                ours = [node((2 * i + int(n[0] + 1), 2 * j + int(n[1] + 1), 2 * k + int(n[2] + 1)))
                        for n in NATURAL]
                gmsh = [0] * 20
                for a, g in enumerate(GMSH_HEX20):
                    gmsh[g] = ours[a]
                elements.append(gmsh)
                if j == 0:
                    bottom.append([ours[a] for a in YMIN])
    out = ["$MeshFormat", "4.1 0 8", "$EndMeshFormat", "$PhysicalNames", "2",
           '2 1 "bottom"', '3 10 "sink"', "$EndPhysicalNames", "$Entities", "0 0 1 1"]
    lo = [min(p[d] for p in nodes) for d in range(3)]
    hi = [max(p[d] for p in nodes) for d in range(3)]
    box = " ".join(f"{v:.9g}" for v in lo + hi)
    out += [f"1 {box} 1 1 0", f"1 {box} 1 10 0", "$EndEntities"]
    n = len(nodes)
    out += ["$Nodes", f"1 {n} 1 {n}", f"3 1 0 {n}"]
    out += [str(t) for t in range(1, n + 1)]
    out += [f"{p[0]:.12g} {p[1]:.12g} {p[2]:.12g}" for p in nodes]
    out += ["$EndNodes"]
    total = len(elements) + len(bottom)
    out += ["$Elements", f"2 {total} 1 {total}", f"2 1 16 {len(bottom)}"]
    tag = 1
    for f in bottom:
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
