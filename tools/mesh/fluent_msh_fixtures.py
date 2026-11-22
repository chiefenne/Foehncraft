"""Writes the Fluent `.msh` test meshes of crates/foehncraft-mesh/tests/data/fluent_msh
from the ASCII OpenFOAM test meshes of crates/foehncraft-mesh/tests/data/openfoam.

Plain Python, no dependencies. Run from the repository root:
    python3 -m tools.mesh.fluent_msh_fixtures

- channel_hex.msh: the hexahedral channel, one section per zone with
  fixed-type quadrilateral faces, the two cell zones as contiguous cell
  ranges, zone names from (45 ...) sections.
- channel_hex_flipped.msh: the same with every face turned around (nodes
  reversed, cells swapped): the reader decides the side by geometry.
- dual_polyhedra.msh: the polyhedral dual mesh, mixed faces (the node count
  leads each face), one cell zone.
"""

import os
import re

from tools.mesh.cgns_ngon_fixture import body, foam_numbers

SRC = "crates/foehncraft-mesh/tests/data/openfoam"
OUT = "crates/foehncraft-mesh/tests/data/fluent_msh"


def read_foam(d):
    pm = f"{SRC}/{d}/constant/polyMesh"
    _, rest = foam_numbers(f"{pm}/points")
    pts = [[float(v) for v in p.split()] for p in re.findall(r"\(([^()]+)\)", rest)]
    _, rest = foam_numbers(f"{pm}/faces")
    faces = [
        [int(v) for v in f.split()] for f in re.findall(r"\d+\s*\(([^()]+)\)", rest)
    ]
    owner = [int(v) for v in foam_numbers(f"{pm}/owner")[1].split()]
    neigh = [int(v) for v in foam_numbers(f"{pm}/neighbour")[1].split()]
    patches = []
    for name, kind, nf, start in re.findall(
        r"(\w+)\s*\{[^}]*?type\s+(\w+);[^}]*?nFaces\s+(\d+);\s*startFace\s+(\d+);",
        body(f"{pm}/boundary"),
    ):
        patches.append((name, kind, int(nf), int(start)))
    zones = []
    zp = f"{pm}/cellZones"
    if os.path.exists(zp):
        for name, n, labels in re.findall(
            r"(\w+)\s*\{[^}]*?cellLabels\s+List<label>\s*(\d+)\s*\(([^)]*)\)", body(zp)
        ):
            zones.append((name, [int(v) for v in labels.split()]))
    return pts, faces, owner, neigh, patches, zones


def write(path, pts, faces, owner, neigh, patches, zones, mixed, flip=False):
    n_cells = max(max(owner), max(neigh)) + 1
    # cells of one zone must be contiguous: number the zones one after the other
    new = list(range(n_cells))
    zone_ranges = []
    if zones:
        order = [c for _, cells in zones for c in cells]
        assert sorted(order) == list(range(n_cells))
        new = [0] * n_cells
        for i, c in enumerate(order):
            new[c] = i
        first = 1
        for name, cells in zones:
            zone_ranges.append((name, first, first + len(cells) - 1))
            first += len(cells)
    else:
        zone_ranges.append(("fluid", 1, n_cells))
    n_inner = len(neigh)
    sections = [("interior-1", "interior", 2, 0, n_inner)]
    for name, kind, nf, start in patches:
        bc = {"wall": 3, "symmetry": 7, "empty": 7}.get(kind, 4)
        label = {"wall": "wall", "symmetry": "symmetry", "empty": "symmetry"}.get(
            kind, "pressure-inlet"
        )
        sections.append((name, label, bc, start, start + nf))

    def f(v):
        return format(v, "x")

    with open(path, "w") as o:
        o.write(
            '(0 "Fluent .msh test mesh, tools/mesh/fluent_msh_fixtures.py")\n(0 "")\n'
        )
        o.write("(2 3)\n")
        o.write(f"(10 (0 1 {f(len(pts))} 0 3))\n")
        o.write(f"(12 (0 1 {f(n_cells)} 0 0))\n")
        o.write(f"(13 (0 1 {f(len(faces))} 0 0))\n")
        o.write(f"(10 (1 1 {f(len(pts))} 1 3)(\n")
        o.writelines("  " + " ".join(f"{v:.16e}" for v in p) + "\n" for p in pts)
        o.write("))\n")
        for i, (name, a, b) in enumerate(zone_ranges):
            o.write(f"(12 ({f(i + 1)} {f(a)} {f(b)} 1 {'7' if mixed else '4'}))\n")
        for z, (name, label, bc, a, b) in enumerate(sections):
            zid = 10 + z
            sizes = {len(faces[i]) for i in range(a, b)}
            ftype = 0 if mixed or len(sizes) != 1 else sizes.pop()
            o.write(f"(13 ({f(zid)} {f(a + 1)} {f(b)} {f(bc)} {f(ftype)})(\n")
            for i in range(a, b):
                nodes = [n + 1 for n in faces[i]]
                c0 = new[owner[i]] + 1
                c1 = new[neigh[i]] + 1 if i < n_inner else 0
                if flip:
                    nodes.reverse()
                    if c1:
                        c0, c1 = c1, c0
                lead = f"{f(len(nodes))} " if ftype == 0 else ""
                o.write(
                    "  " + lead + " ".join(f(n) for n in nodes) + f" {f(c0)} {f(c1)}\n"
                )
            o.write("))\n")
        for i, (name, a, b) in enumerate(zone_ranges):
            o.write(f"(45 ({f(i + 1)} fluid {name})())\n")
        for z, (name, label, bc, a, b) in enumerate(sections):
            o.write(f"(45 ({f(10 + z)} {label} {name})())\n")


def main():
    os.makedirs(OUT, exist_ok=True)
    pts, faces, owner, neigh, patches, zones = read_foam("channel_ascii")
    write(
        f"{OUT}/channel_hex.msh", pts, faces, owner, neigh, patches, zones, mixed=False
    )
    write(
        f"{OUT}/channel_hex_flipped.msh",
        pts,
        faces,
        owner,
        neigh,
        patches,
        zones,
        mixed=False,
        flip=True,
    )
    pts, faces, owner, neigh, patches, zones = read_foam("dual_polyhedra")
    write(
        f"{OUT}/dual_polyhedra.msh",
        pts,
        faces,
        owner,
        neigh,
        patches,
        zones,
        mixed=True,
    )


if __name__ == "__main__":
    main()
