"""Pictures of the reed valve beam of the structural solver (gallery entry
of MS5, L53).

reed_modes.png: the first four modes of the rigidly clamped reed
(examples/cases/structure/reed_beam.toml) drawn as the reed's 6 mm wide
strip: every node of the beam carries its section, moved by the
translation of the axis and turned by its rotation, so flapping and
torsion both show. Colour: the size of the displacement of the strip
(scaled mode shape), Rainbow Desaturated.

reed_clamp.png: the first two flapping frequencies of the reed against the
rotational clamp spring (D290), from a soft mounting to a rigid one: the
solver (marks, `structure.run`) on the frequency equation of the beam with
a rotational root spring (line), and the rigid clamp (dashed).

    pip install matplotlib numpy
    foehn call structure.run --input examples/cases/tools/structure_reed_beam.json
    python3 -m tools.plot.structure_reed modes output/reed_beam docs/figures/gallery/structure_reed_modes.png
    python3 -m tools.plot.structure_reed clamp target/ci/foehn docs/figures/gallery/structure_reed_clamp.png
"""
import json
import os
import subprocess
import sys
import tempfile

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from tools.plot.colormaps import colormap

RAINBOW = colormap()

CASE = "examples/cases/structure/reed_beam.toml"
E, NU, RHO = 206e9, 0.3, 7850.0
L, B, H = 0.03, 6e-3, 0.2e-3


def read_vtk(path):
    """Points, lines and the point vectors of a legacy VTK polydata file."""
    with open(path) as f:
        tok = f.read().split("\n")
    pts, lines, vec = None, [], {}
    i = 0
    while i < len(tok):
        t = tok[i].split()
        if t and t[0] == "POINTS":
            n = int(t[1])
            pts = np.array([[float(v) for v in tok[i + 1 + k].split()] for k in range(n)])
            i += n
        elif t and t[0] == "LINES":
            m = int(t[1])
            lines = [[int(v) for v in tok[i + 1 + k].split()[1:]] for k in range(m)]
            i += m
        elif t and t[0] == "VECTORS":
            n = len(pts)
            vec[t[1]] = np.array([[float(v) for v in tok[i + 1 + k].split()] for k in range(n)])
            i += n
        i += 1
    return pts, lines, vec


def modes(out_dir, fig_path):
    pts, _, vec = read_vtk(os.path.join(out_dir, "structure.vtk"))
    with open(os.path.join(out_dir, "result.json")) as f:
        res = json.load(f)
    want = [f"mode_{k}" for k in (1, 2, 3, 4)]
    missing = [w for w in want if w not in vec]
    if missing:
        sys.exit(f"{out_dir}/structure.vtk lacks {missing}; it has: {', '.join(sorted(vec))}")
    order = np.argsort(pts[:, 0])
    x = pts[order]
    fig = plt.figure(figsize=(9.0, 6.6))
    for k, name in enumerate(want):
        u = vec[name][order]
        r = vec[name + "_rotation"][order]
        # the strip: section points at y = ±B/2 and in between, moved by
        # u + θ × (0, y, 0), with the displacement scaled to 0.35 L
        # (a torsion mode moves the strip's edges, not its axis: the scale
        # is that of the largest displacement of the strip)
        ys = np.linspace(-B / 2, B / 2, 7)
        D = np.stack([u + np.cross(r, np.array([0.0, y, 0.0])) for y in ys], axis=1)
        C = np.linalg.norm(D, axis=2)
        big = C.max()
        C /= big
        X = np.zeros((len(x), len(ys), 3))
        for j, y in enumerate(ys):
            X[:, j, :] = x + np.array([0.0, y, 0.0]) + 0.3 * L * D[:, j, :] / big
        ax = fig.add_subplot(2, 2, k + 1, projection="3d")
        quads, cols = [], []
        for a in range(len(x) - 1):
            for j in range(len(ys) - 1):
                quads.append([X[a, j], X[a + 1, j], X[a + 1, j + 1], X[a, j + 1]])
                cols.append(0.25 * (C[a, j] + C[a + 1, j] + C[a + 1, j + 1] + C[a, j + 1]))
        pc = Poly3DCollection(quads, facecolors=RAINBOW(np.array(cols)), edgecolors="none")
        ax.add_collection3d(pc)
        # the undeformed reed and the clamp
        ax.plot([0, L, L, 0, 0], [-B / 2, -B / 2, B / 2, B / 2, -B / 2], [0, 0, 0, 0, 0],
                color="0.55", lw=0.8, ls=":")
        ax.plot([0, 0], [-B / 2 - 1e-3, B / 2 + 1e-3], [0, 0], color="k", lw=3)
        mo = res["modal"]["modes"][k]
        kind = "torsion" if np.abs(u[:, 2]).max() < 0.2 * np.abs(np.cross(r, [0, B / 2, 0])).max() else "flapping"
        ax.set_title(f"mode {k + 1}: {mo['frequency']:.0f} Hz ({kind})", fontsize=10)
        ax.set_xlim(0, L)
        ax.set_ylim(-0.35 * L, 0.35 * L)
        ax.set_zlim(-0.35 * L, 0.35 * L)
        ax.set_box_aspect((1.0, 0.7, 0.7))
        ax.view_init(elev=22, azim=-58)
        ax.set_xticks([0, 0.01, 0.02, 0.03])
        ax.set_xticklabels(["0", "10", "20", "30 mm"], fontsize=7)
        ax.set_yticks([])
        ax.set_zticks([])
        ax.xaxis.pane.set_alpha(0.0)
        ax.yaxis.pane.set_alpha(0.0)
        ax.zaxis.pane.set_alpha(0.0)
    sm = plt.cm.ScalarMappable(cmap=RAINBOW)
    cb = fig.colorbar(sm, ax=fig.axes, shrink=0.6, pad=0.02)
    cb.set_label("displacement / largest displacement of the mode")
    fig.suptitle("Reed valve, 30 × 6 × 0.2 mm spring steel, rigid clamp: the first four modes", fontsize=11)
    fig.savefig(fig_path, dpi=110)
    print(fig_path)


def ei():
    return E * B * H**3 / 12.0


def spring_roots(k, count=2):
    """The first `count` roots βL of the beam with a rotational root
    spring k (N·m/rad), free at the tip."""

    def det(bl):
        beta = bl / L
        ch, sh, c, s = np.cosh(bl), np.sinh(bl), np.cos(bl), np.sin(bl)
        m = np.array([[2 * ei() * beta, -k, -k], [ch + c, sh, -s], [sh - s, ch, -c]])
        return np.linalg.det(m)

    roots, a = [], 1e-3
    fa = det(a)
    while len(roots) < count:
        b = a + 1e-3
        fb = det(b)
        if np.sign(fa) != np.sign(fb):
            lo, hi, flo = a, b, fa
            for _ in range(80):
                mid = 0.5 * (lo + hi)
                fm = det(mid)
                if np.sign(fm) == np.sign(flo):
                    lo, flo = mid, fm
                else:
                    hi = mid
            roots.append(0.5 * (lo + hi))
        a, fa = b, fb
    return roots


def freq(bl):
    return bl**2 / L**2 * np.sqrt(ei() / (RHO * B * H)) / (2 * np.pi)


def flapping(res, shapes):
    """Frequencies of the modes that are mostly w (not torsion or in-plane)."""
    out = []
    for mo, (u, r) in zip(res["modal"]["modes"], shapes):
        w = np.sum(u[:, 2] ** 2)
        if w > np.sum(u[:, 1] ** 2) and w > np.sum((r[:, 0] * B) ** 2):
            out.append(mo["frequency"])
    return out


def clamp(foehn, fig_path):
    with open(CASE) as f:
        base = f.read().replace('theory = "timoshenko"', 'theory = "euler-bernoulli"')
    base = base.replace("static = true", "static = false")
    ratios = np.logspace(-1, 3, 13)  # k L / EI
    fe = []
    with tempfile.TemporaryDirectory() as tmp:
        for q in ratios:
            k = q * ei() / L / B  # per metre of clamp line
            text = base.replace('line = "y"', f'line = "y"\nk_rotation = {float(k)!r}')
            case = os.path.join(tmp, "reed.toml")
            with open(case, "w") as f:
                f.write(text)
            inp = os.path.join(tmp, "in.json")
            out = os.path.join(tmp, "out")
            with open(inp, "w") as f:
                json.dump({"case": case, "out": out}, f)
            subprocess.run([foehn, "call", "structure.run", "--input", inp, "--quiet"],
                           check=True, stdout=subprocess.DEVNULL)
            with open(os.path.join(out, "result.json")) as f:
                res = json.load(f)
            _, _, vec = read_vtk(os.path.join(out, "structure.vtk"))
            n = len(res["modal"]["modes"])
            shapes = [(vec[f"mode_{i}"], vec[f"mode_{i}_rotation"]) for i in range(1, n + 1)]
            fe.append(flapping(res, shapes)[:2])
    fe = np.array(fe)
    dense = np.logspace(-1, 3, 200)
    exact = np.array([[freq(b) for b in spring_roots(q * ei() / L)] for q in dense])
    rigid = [freq(1.875104068711961), freq(4.694091132974175)]
    fig, ax = plt.subplots(figsize=(7.4, 4.6))
    colours = [RAINBOW(0.12), RAINBOW(0.85)]
    for m in range(2):
        ax.plot(dense, exact[:, m], color=colours[m], lw=1.6,
                label=f"mode {m + 1}, frequency equation")
        ax.plot(ratios, fe[:, m], "o", ms=5, mfc="white", mec=colours[m],
                label=f"mode {m + 1}, structure.run (40 beams)")
        ax.axhline(rigid[m], color=colours[m], ls="--", lw=0.9)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("rotational clamp stiffness  k L / EI")
    ax.set_ylabel("natural frequency (Hz)")
    ax.set_title("Reed on a rotational clamp spring: flapping modes (dashed: rigid clamp)", fontsize=10)
    ax.grid(True, which="both", lw=0.3, color="0.85")
    ax.legend(fontsize=8, loc="lower right")
    err = np.max(np.abs(fe - np.array([[freq(b) for b in spring_roots(q * ei() / L)] for q in ratios])) / fe)
    print(f"largest relative difference to the frequency equation: {err:.2e}")
    fig.tight_layout()
    fig.savefig(fig_path, dpi=110)
    print(fig_path)


if __name__ == "__main__":
    if len(sys.argv) != 4 or sys.argv[1] not in ("modes", "clamp"):
        sys.exit(__doc__)
    {"modes": modes, "clamp": clamp}[sys.argv[1]](sys.argv[2], sys.argv[3])
