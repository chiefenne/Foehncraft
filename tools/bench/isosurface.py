"""Iso-surface extraction benchmark (L13 stage 1, D226; owner rule: benchmark
the computer-graphics methods).

Marching tetrahedra (Doi & Koide 1991, what `foehncraft_mesh::isosurface`
implements) against dual contouring (Ju, Losasso, Schaefer & Warren 2002)
and subgrid marching tetrahedra (Baktash, Gillespie & Crane 2026, read in
full: ACM TOG 45(4) 57, CC BY 4.0) on fields with an exact surface:

    sphere      signed distance of a sphere (smooth)
    box         signed distance of a rotated box (sharp edges and corners)
    plate       signed distance of a plate thinner than a cell (sub-grid)
    fraction    a smeared solid fraction of the sphere, clamp(½ + d/h): what
                a penalised body or a design density gives

on a grid of n³ cubes (each split into 6 tetrahedra for marching tets). The
metrics: two-sided Hausdorff distance to the exact surface (in cells),
area and volume error, edges that are open or non-manifold, triangle pairs
that intersect, the smallest triangle angle and the time.

Subgrid marching needs every crossing of the surface along each grid edge.
For data given at the nodes and interpolated linearly along the edges (the
cell and node values of a CFD mesh) an edge has at most one crossing, so
subgrid marching gives exactly the marching-tets surface. Its gain appears
when the field itself can be evaluated along an edge (an analytic signed
distance, stage 0 and stage 3): the script counts the grid edges with two or
more crossings of the exact field (sampled at 64 points per edge), which
marching tets and dual contouring miss and subgrid marching keeps.

    pip install numpy matplotlib
    python3 -m tools.bench.isosurface               the table (markdown)
    python3 -m tools.bench.isosurface --json OUT    the numbers as JSON
    python3 -m tools.bench.isosurface --figure PNG  the gallery image
"""
import argparse
import json
import sys
import time

import numpy as np

# the six tetrahedra of a cube around its diagonal 0–7 (corners by bits x, y, z)
KUHN = [(0, 1, 3, 7), (0, 3, 2, 7), (0, 2, 6, 7), (0, 6, 4, 7), (0, 4, 5, 7), (0, 5, 1, 7)]
CORNERS = np.array([[(c >> 0) & 1, (c >> 1) & 1, (c >> 2) & 1] for c in range(8)])


def rotation(ax, az):
    cx, sx, cz, sz = np.cos(ax), np.sin(ax), np.cos(az), np.sin(az)
    rx = np.array([[1, 0, 0], [0, cx, -sx], [0, sx, cx]])
    rz = np.array([[cz, -sz, 0], [sz, cz, 0], [0, 0, 1]])
    return rz @ rx


def sdf_box(x, half, rot=None):
    """Signed distance of a box (negative inside), exact."""
    p = x if rot is None else x @ rot  # world -> box frame (rot orthogonal)
    q = np.abs(p) - half
    out = np.linalg.norm(np.maximum(q, 0.0), axis=-1)
    return out + np.minimum(np.max(q, axis=-1), 0.0)


def shapes(h):
    rot = rotation(np.radians(20.0), np.radians(30.0))
    r = 0.6
    box_half = np.array([0.45, 0.45, 0.45])
    plate_half = np.array([0.7, 0.7, 0.15 * h])
    sphere = lambda x: np.linalg.norm(x, axis=-1) - r  # noqa: E731
    return {
        "sphere": dict(
            f=sphere, d=sphere, area=4 * np.pi * r * r, volume=4 / 3 * np.pi * r**3,
            sample=lambda n: sample_sphere(r, n),
        ),
        "box": dict(
            f=lambda x: sdf_box(x, box_half, rot), d=lambda x: sdf_box(x, box_half, rot),
            area=6 * 0.9**2, volume=0.9**3,
            sample=lambda n: sample_box(box_half, rot, n),
        ),
        "plate": dict(
            f=lambda x: sdf_box(x - [0, 0, 0.013], plate_half),
            d=lambda x: sdf_box(x - [0, 0, 0.013], plate_half),
            area=2 * 1.4 * 1.4 + 4 * 1.4 * 0.3 * h, volume=1.4 * 1.4 * 0.3 * h,
            sample=lambda n: sample_box(plate_half, None, n) + [0, 0, 0.013],
        ),
        # a smeared solid fraction (positive = solid here: level ½ of
        # 1 − fraction is the same surface), its exact surface the sphere
        "fraction": dict(
            f=lambda x: -(np.clip(0.5 - sphere(x) / h, 0.0, 1.0) - 0.5), d=sphere,
            area=4 * np.pi * r * r, volume=4 / 3 * np.pi * r**3,
            sample=lambda n: sample_sphere(r, n),
        ),
    }


def sample_sphere(r, n):
    k = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * k / n)
    th = np.pi * (1 + 5**0.5) * k
    return r * np.stack([np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)], -1)


def sample_box(half, rot, n):
    m = max(int(np.sqrt(n / 6)), 2)
    u = (np.arange(m) + 0.5) / m * 2 - 1
    a, b = np.meshgrid(u, u)
    a, b = a.ravel(), b.ravel()
    pts = []
    for ax in range(3):
        o = [i for i in range(3) if i != ax]
        for s in (-1, 1):
            p = np.zeros((a.size, 3))
            p[:, ax] = s * half[ax]
            p[:, o[0]] = a * half[o[0]]
            p[:, o[1]] = b * half[o[1]]
            pts.append(p)
    p = np.concatenate(pts)
    return p if rot is None else p @ rot.T


def grid(n):
    x = np.linspace(-1, 1, n + 1)
    return np.stack(np.meshgrid(x, x, x, indexing="ij"), -1), 2.0 / n


# --------------------------------------------------------------------------
# marching tetrahedra


def marching_tets(pts, val):
    n = val.shape[0] - 1
    idx = np.arange((n + 1) ** 3).reshape((n + 1,) * 3)
    base = idx[:-1, :-1, :-1].ravel()
    off = np.array([idx[c[0], c[1], c[2]] for c in CORNERS]) - idx[0, 0, 0]
    P = pts.reshape(-1, 3)
    V = val.ravel()
    verts, key_of, tris = [], {}, []

    def vertex(a, b):
        k = (min(a, b), max(a, b))
        if k not in key_of:
            p, q = k
            t = V[p] / (V[p] - V[q])
            key_of[k] = len(verts)
            verts.append(P[p] + t * (P[q] - P[p]))
        return key_of[k]

    corner = base[:, None] + off[None, :]
    for t in KUHN:
        g = corner[:, list(t)]  # (cubes, 4)
        inside = V[g] < 0  # negative inside
        cnt = inside.sum(1)
        for c in np.nonzero((cnt > 0) & (cnt < 4))[0]:
            ids = g[c]
            ins = [ids[k] for k in range(4) if inside[c, k]]
            out = [ids[k] for k in range(4) if not inside[c, k]]
            if len(ins) == 1:
                poly = [vertex(ins[0], o) for o in out]
            elif len(ins) == 3:
                poly = [vertex(i, out[0]) for i in ins]
            else:
                poly = [vertex(ins[0], out[0]), vertex(ins[0], out[1]),
                        vertex(ins[1], out[1]), vertex(ins[1], out[0])]
            d = P[out].mean(0) - P[ins].mean(0)
            pv = np.array([verts[k] for k in poly])
            nrm = sum(np.cross(pv[k] - pv[0], pv[k + 1] - pv[0]) for k in range(1, len(pv) - 1))
            flip = np.dot(nrm, d) < 0
            for k in range(1, len(poly) - 1):
                tri = [poly[0], poly[k + 1], poly[k]] if flip else [poly[0], poly[k], poly[k + 1]]
                if len(set(tri)) == 3:
                    tris.append(tri)
    if not tris:
        return np.zeros((0, 3)), np.zeros((0, 3), dtype=np.int64)
    return np.array(verts), np.array(tris, dtype=np.int64)


# --------------------------------------------------------------------------
# dual contouring


def gradient(f, x, h):
    e = np.eye(3) * (1e-3 * h)
    g = np.stack([(f(x + e[i]) - f(x - e[i])) / (2e-3 * h) for i in range(3)], -1)
    return g / np.maximum(np.linalg.norm(g, axis=-1, keepdims=True), 1e-300)


def dual_contouring(pts, val, f, h):
    """Dual contouring: one vertex per sign-changing cube at the minimiser of
    the QEF of its edge crossings and normals (pulled to the mass point by a
    small regularisation, clamped to the cube); one quad per crossed edge."""
    n = val.shape[0] - 1
    vid = -np.ones((n, n, n), dtype=np.int64)
    verts = []
    edges = []
    for ax in range(3):
        sl0 = [slice(None)] * 3
        sl1 = [slice(None)] * 3
        sl0[ax] = slice(0, -1)
        sl1[ax] = slice(1, None)
        a, b = val[tuple(sl0)], val[tuple(sl1)]
        cross = (a < 0) != (b < 0)
        for i, j, k in np.argwhere(cross):
            p0 = np.array([i, j, k])
            p1 = p0.copy()
            p1[ax] += 1
            va, vb = val[tuple(p0)], val[tuple(p1)]
            t = va / (va - vb)
            x = pts[tuple(p0)] + t * (pts[tuple(p1)] - pts[tuple(p0)])
            edges.append((ax, p0, x, va < 0))
    if not edges:
        return np.zeros((0, 3)), np.zeros((0, 3), dtype=np.int64)
    xs = np.array([e[2] for e in edges])
    ns = gradient(f, xs, h)
    per_cube = {}
    for (ax, p0, x, _), nn in zip(edges, ns):
        others = [d for d in range(3) if d != ax]
        for da in (0, -1):
            for db in (0, -1):
                c = p0.copy()
                c[others[0]] += da
                c[others[1]] += db
                if np.all(c >= 0) and np.all(c < n):
                    per_cube.setdefault(tuple(c), []).append((x, nn))
    lo = pts[0, 0, 0]
    for c in sorted(per_cube):
        data = per_cube[c]
        X = np.array([d[0] for d in data])
        N = np.array([d[1] for d in data])
        m = X.mean(0)
        lam = 0.05
        A = np.vstack([N, np.sqrt(lam) * np.eye(3)])
        rhs = np.concatenate([np.sum(N * (X - m), 1), np.zeros(3)])
        y = np.linalg.lstsq(A, rhs, rcond=None)[0] + m
        cmin = lo + np.array(c) * h
        y = np.clip(y, cmin, cmin + h)
        vid[c] = len(verts)
        verts.append(y)
    tris = []
    for ax, p0, x, neg_first in edges:
        others = [d for d in range(3) if d != ax]
        quad = []
        for da, db in ((0, 0), (-1, 0), (-1, -1), (0, -1)):
            c = p0.copy()
            c[others[0]] += da
            c[others[1]] += db
            if np.any(c < 0) or np.any(c >= n):
                quad = None
                break
            quad.append(vid[tuple(c)])
        if quad is None:
            continue
        # orientation: the normal out of the solid (negative inside), that
        # is along +ax when the edge starts inside
        q = np.array([verts[k] for k in quad])
        nrm = np.cross(q[1] - q[0], q[2] - q[0]) + np.cross(q[2] - q[0], q[3] - q[0])
        if (nrm[ax] > 0) != neg_first:
            quad = quad[::-1]
        tris.append([quad[0], quad[1], quad[2]])
        tris.append([quad[0], quad[2], quad[3]])
    return np.array(verts), np.array(tris, dtype=np.int64)


# --------------------------------------------------------------------------
# metrics


def area_volume(V, T):
    a, b, c = V[T[:, 0]], V[T[:, 1]], V[T[:, 2]]
    cr = np.cross(b - a, c - a)
    return 0.5 * np.linalg.norm(cr, axis=1).sum(), np.einsum("ij,ij->i", a, np.cross(b, c)).sum() / 6


def edge_counts(T):
    e = np.sort(np.concatenate([T[:, [0, 1]], T[:, [1, 2]], T[:, [2, 0]]]), axis=1)
    _, cnt = np.unique(e, axis=0, return_counts=True)
    return int((cnt == 1).sum()), int((cnt > 2).sum())


def min_angles(V, T):
    p = V[T]
    out = np.full(len(T), 180.0)
    for k in range(3):
        a = p[:, (k + 1) % 3] - p[:, k]
        b = p[:, (k + 2) % 3] - p[:, k]
        c = np.einsum("ij,ij->i", a, b) / (np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1))
        out = np.minimum(out, np.degrees(np.arccos(np.clip(c, -1, 1))))
    return out


def point_triangle_distance(P, A, B, C):
    """Distance from each point to each triangle (Ericson 2005, 5.1.5),
    vectorised over the broadcast shape."""
    ab, ac, ap = B - A, C - A, P - A
    d1, d2 = np.sum(ab * ap, -1), np.sum(ac * ap, -1)
    bp = P - B
    d3, d4 = np.sum(ab * bp, -1), np.sum(ac * bp, -1)
    cp = P - C
    d5, d6 = np.sum(ab * cp, -1), np.sum(ac * cp, -1)
    va = d3 * d6 - d5 * d4
    vb = d5 * d2 - d1 * d6
    vc = d1 * d4 - d3 * d2
    denom = np.where(va + vb + vc == 0, 1e-300, va + vb + vc)
    v = vb / denom
    w = vc / denom
    Q = A + ab * v[..., None] + ac * w[..., None]
    # regions
    def put(mask, X):
        nonlocal Q
        Q = np.where(mask[..., None], X, Q)
    with np.errstate(divide="ignore", invalid="ignore"):
        e_bc = (va <= 0) & (d4 - d3 >= 0) & (d5 - d6 >= 0)
        put(e_bc, B + (C - B) * ((d4 - d3) / ((d4 - d3) + (d5 - d6)))[..., None])
        e_ac = (vb <= 0) & (d2 >= 0) & (d6 <= 0)
        put(e_ac, A + ac * (d2 / (d2 - d6))[..., None])
        e_ab = (vc <= 0) & (d1 >= 0) & (d3 <= 0)
        put(e_ab, A + ab * (d1 / (d1 - d3))[..., None])
        put((d6 >= 0) & (d5 <= d6), C)
        put((d3 >= 0) & (d4 <= d3), B)
        put((d1 <= 0) & (d2 <= 0), A)
    return np.linalg.norm(P - Q, axis=-1)


def hausdorff(V, T, d, sample):
    # mesh -> surface: |d| at points on the triangles (vertices, midpoints, centroid)
    a, b, c = V[T[:, 0]], V[T[:, 1]], V[T[:, 2]]
    pts = np.concatenate([V, (a + b) / 2, (b + c) / 2, (c + a) / 2, (a + b + c) / 3])
    m2s = np.abs(d(pts)).max()
    # surface -> mesh: points of the exact surface against all triangles
    S = sample(1500)
    best = np.full(len(S), np.inf)
    for k in range(0, len(T), 800):
        A, B, C = (x[None, k:k + 800] for x in (a, b, c))
        best = np.minimum(best, point_triangle_distance(S[:, None], A, B, C).min(1))
    return max(m2s, best.max())


def segment_triangle(P0, P1, A, B, C, eps=1e-12):
    """True where the segment P0–P1 crosses the triangle ABC (strictly)."""
    e1, e2 = B - A, C - A
    d = P1 - P0
    h = np.cross(d, e2)
    det = np.sum(e1 * h, -1)
    ok = np.abs(det) > eps
    inv = np.where(ok, 1 / np.where(ok, det, 1), 0)
    s = P0 - A
    u = inv * np.sum(s * h, -1)
    q = np.cross(s, e1)
    v = inv * np.sum(d * q, -1)
    t = inv * np.sum(e2 * q, -1)
    tol = 1e-9
    return ok & (u > tol) & (v > tol) & (u + v < 1 - tol) & (t > tol) & (t < 1 - tol)


def self_intersections(V, T):
    """Pairs of triangles without a common vertex that intersect (an edge of
    one crosses the other), candidates by overlapping bounding boxes."""
    lo, hi = V[T].min(1), V[T].max(1)
    order = np.argsort(lo[:, 0])
    lo, hi, Ts = lo[order], hi[order], T[order]
    count = 0
    n = len(Ts)
    for i in range(n):
        j0 = i + 1
        j1 = np.searchsorted(lo[:, 0], hi[i, 0], side="right")
        if j1 <= j0:
            continue
        js = np.arange(j0, j1)
        js = js[np.all((lo[js] <= hi[i]) & (hi[js] >= lo[i]), axis=1)]
        if js.size == 0:
            continue
        share = np.any(np.isin(Ts[js], Ts[i]), axis=1)
        js = js[~share]
        if js.size == 0:
            continue
        a = V[Ts[i]]
        b = V[Ts[js]]
        hit = np.zeros(js.size, bool)
        for k in range(3):
            hit |= segment_triangle(a[k], a[(k + 1) % 3], b[:, 0], b[:, 1], b[:, 2])
            hit |= segment_triangle(b[:, k], b[:, (k + 1) % 3], a[0], a[1], a[2])
        count += int(hit.sum())
    return count


def multi_crossings(pts, f, n, samples=64):
    """Grid edges along which the exact field changes sign twice or more."""
    total = 0
    s = np.linspace(0, 1, samples + 1)
    for ax in range(3):
        sl0 = [slice(None)] * 3
        sl0[ax] = slice(0, -1)
        p0 = pts[tuple(sl0)].reshape(-1, 3)
        step = np.zeros(3)
        step[ax] = 2.0 / n
        for k in range(0, len(p0), 4096):
            q = p0[k:k + 4096, None, :] + s[None, :, None] * step
            v = f(q)
            changes = np.sum(np.signbit(v[:, 1:]) != np.signbit(v[:, :-1]), axis=1)
            total += int((changes >= 2).sum())
    return total


def run(n=24):
    pts, h = grid(n)
    rows = []
    surfaces = {}
    for name, sh in shapes(h).items():
        val = sh["f"](pts)
        for method in ("marching tets", "dual contouring"):
            t0 = time.perf_counter()
            if method == "marching tets":
                V, T = marching_tets(pts, val)
            else:
                V, T = dual_contouring(pts, val, sh["f"], h)
            dt = time.perf_counter() - t0
            surfaces[(name, method)] = (V, T)
            if len(T) == 0 or len(V) == 0:
                rows.append(dict(field=name, method=method, triangles=0, missed=True))
                continue
            area, vol = area_volume(V, T)
            open_e, nonman = edge_counts(T)
            ang = min_angles(V, T)
            rows.append(dict(
                field=name, method=method, triangles=int(len(T)),
                hausdorff_cells=float(hausdorff(V, T, sh["d"], sh["sample"]) / h),
                area_error=float(area / sh["area"] - 1), volume_error=float(vol / sh["volume"] - 1),
                open_edges=open_e, nonmanifold_edges=nonman,
                intersections=self_intersections(V, T),
                min_angle=float(ang.min()), slivers=int((ang < 5).sum()), time_s=dt,
            ))
        if name != "fraction":
            rows.append(dict(field=name, method="subgrid marching tets (edges with ≥ 2 crossings)",
                             multi_crossing_edges=multi_crossings(pts, sh["f"], n)))
    return dict(n=n, h=h, rows=rows), surfaces


def table(res):
    out = [f"Grid {res['n']}³ cubes on [−1, 1]³ (h = {res['h']:.4f}); Hausdorff in cells.", "",
           "| field | method | triangles | Hausdorff | area | volume | open | non-manifold | "
           "intersecting pairs | min angle | < 5° | time (s) |",
           "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in res["rows"]:
        if "multi_crossing_edges" in r:
            out.append(f"| {r['field']} | {r['method']} | grid edges: {r['multi_crossing_edges']} "
                       "| | | | | | | | | |")
        elif r.get("missed"):
            out.append(f"| {r['field']} | {r['method']} | 0 (surface missed) | | | | | | | | | |")
        else:
            out.append(
                f"| {r['field']} | {r['method']} | {r['triangles']} | {r['hausdorff_cells']:.3f} | "
                f"{100 * r['area_error']:+.2f} % | {100 * r['volume_error']:+.2f} % | {r['open_edges']} | "
                f"{r['nonmanifold_edges']} | {r['intersections']} | {r['min_angle']:.1f}° | "
                f"{r['slivers']} | {r['time_s']:.2f} |")
    return "\n".join(out)


def figure(surfaces, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    fig = plt.figure(figsize=(12, 6.4), dpi=150)
    panels = [("box", "marching tets"), ("box", "dual contouring"),
              ("fraction", "marching tets"), ("plate", "marching tets")]
    titles = ["rotated box: marching tetrahedra", "rotated box: dual contouring",
              "smeared solid fraction of a sphere:\nmarching tetrahedra",
              "plate 0.3 cells thick:\nmarching tetrahedra miss it"]
    light = np.array([0.4, -0.5, 0.75])
    light /= np.linalg.norm(light)
    cmap = plt.get_cmap("viridis")
    for k, ((name, method), title) in enumerate(zip(panels, titles)):
        ax = fig.add_subplot(1, 4, k + 1, projection="3d")
        V, T = surfaces[(name, method)]
        if len(T):
            tri = V[T]
            nrm = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
            nrm /= np.maximum(np.linalg.norm(nrm, axis=1, keepdims=True), 1e-300)
            shade = 0.25 + 0.75 * np.abs(nrm @ light)
            col = cmap(0.15 + 0.7 * (tri[:, :, 2].mean(1) + 1) / 2)
            col[:, :3] *= shade[:, None]
            ax.add_collection3d(Poly3DCollection(tri, facecolors=col, edgecolors=(0, 0, 0, 0.25),
                                                 linewidths=0.15))
        else:
            ax.text(0, 0, 0, "no triangles", ha="center")
        for s in "xyz":
            getattr(ax, f"set_{s}lim")(-0.8, 0.8)
        ax.set_box_aspect((1, 1, 1))
        ax.view_init(elev=24, azim=-58)
        ax.set_axis_off()
        ax.set_title(title, fontsize=10)
    fig.suptitle("Iso-surface extraction benchmark (24³ grid): marching tetrahedra, dual contouring",
                 fontsize=12)
    fig.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=0.86, wspace=0.0)
    fig.savefig(path, facecolor="white")
    print(f"wrote {path}")


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--n", type=int, default=24)
    ap.add_argument("--json")
    ap.add_argument("--figure")
    a = ap.parse_args(argv)
    res, surfaces = run(a.n)
    print(table(res))
    if a.json:
        with open(a.json, "w") as f:
            json.dump(res, f, indent=1)
    if a.figure:
        figure(surfaces, a.figure)


if __name__ == "__main__":
    main(sys.argv[1:])
