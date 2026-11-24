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
- channel_hex_binary.msh: channel_hex.msh with binary single-precision
  sections (2010, 2013): little-endian 32-bit reals and integers.
- dual_polyhedra_binary.msh: dual_polyhedra.msh with binary double-precision
  sections (3010, 3012 with the element type of every cell, 3013 mixed
  faces) and a skipped binary section the reader does not need (3059).
"""

import os
import re
import struct

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


def binary_section(o, index, header, fmt, values):
    """A binary section: text header, little-endian data, closing words."""
    o.write(f"({index} ({header})(".encode())
    o.write(struct.pack(f"<{len(values)}{fmt}", *values))
    o.write(f")\nEnd of Binary Section {index})\n".encode())


def write(
    path, pts, faces, owner, neigh, patches, zones, mixed, flip=False, binary=None
):
    """`binary`: None (ASCII), "single" (20xx sections) or "double" (30xx)."""
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

    if binary:
        prefix = 2000 if binary == "single" else 3000
        with open(path, "wb") as o:
            write_binary(
                o, prefix, pts, faces, owner, neigh, new, zone_ranges, sections, mixed
            )
        return
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


def write_binary(
    o, prefix, pts, faces, owner, neigh, new, zone_ranges, sections, mixed
):
    def f(v):
        return format(v, "x")

    n_cells = len(new)
    n_inner = len(neigh)
    o.write(b'(0 "Fluent .msh binary test mesh, tools/mesh/fluent_msh_fixtures.py")\n')
    o.write(b"(2 3)\n")
    o.write(f"(10 (0 1 {f(len(pts))} 0 3))\n".encode())
    o.write(f"(12 (0 1 {f(n_cells)} 0 0))\n".encode())
    o.write(f"(13 (0 1 {f(len(faces))} 0 0))\n".encode())
    real = "f" if prefix == 2000 else "d"
    binary_section(
        o,
        prefix + 10,
        f"1 1 {f(len(pts))} 1 3",
        real,
        [v for p in pts for v in p],
    )
    for i, (name, a, b) in enumerate(zone_ranges):
        if mixed:
            # mixed element type 0: the type of every cell follows (7 = polyhedron)
            binary_section(
                o, prefix + 12, f"{f(i + 1)} {f(a)} {f(b)} 1 0", "i", [7] * (b - a + 1)
            )
        else:
            o.write(f"({prefix + 12} ({f(i + 1)} {f(a)} {f(b)} 1 4))\n".encode())
    for z, (name, label, bc, a, b) in enumerate(sections):
        zid = 10 + z
        sizes = {len(faces[i]) for i in range(a, b)}
        ftype = 0 if mixed or len(sizes) != 1 else sizes.pop()
        ints = []
        for i in range(a, b):
            nodes = [n + 1 for n in faces[i]]
            c0 = new[owner[i]] + 1
            c1 = new[neigh[i]] + 1 if i < n_inner else 0
            if ftype == 0:
                ints.append(len(nodes))
            ints += nodes + [c0, c1]
        binary_section(
            o, prefix + 13, f"{f(zid)} {f(a + 1)} {f(b)} {f(bc)} {f(ftype)}", "i", ints
        )
    if mixed:
        # a section the reader skips; its data hold ')' and '(' on purpose
        binary_section(o, prefix + 59, "1 1 2 1", "i", [0x29, 0x28, 41, 40])
    for i, (name, a, b) in enumerate(zone_ranges):
        o.write(f"(45 ({f(i + 1)} fluid {name})())\n".encode())
    for z, (name, label, bc, a, b) in enumerate(sections):
        o.write(f"(45 ({f(10 + z)} {label} {name})())\n".encode())


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
    write(
        f"{OUT}/channel_hex_binary.msh",
        pts,
        faces,
        owner,
        neigh,
        patches,
        zones,
        mixed=False,
        binary="single",
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
    write(
        f"{OUT}/dual_polyhedra_binary.msh",
        pts,
        faces,
        owner,
        neigh,
        patches,
        zones,
        mixed=True,
        binary="double",
    )


if __name__ == "__main__":
    main()
