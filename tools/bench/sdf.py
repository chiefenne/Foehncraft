"""Signed distance benchmark (L13 stage 0, D225 item 4; decision Q11: the
default method is decided by this benchmark; owner rule: benchmark the
computer-graphics methods).

Reads the fields of `examples/sdf_benchmark.rs` (ours: exact closest point
+ winding-number sign; sweep: fast sweeping + scan-conversion sign; heat:
the signed heat method of Feng & Crane 2024) and compares them with the
exact signed distance of the shape (the clean sphere for the dirty
variants):

    error    largest and RMS |φ − φ_exact| over the nodes within 0.1 of the
             surface (near) and over all nodes, in cells
    sign     nodes farther than one cell from the surface whose sign is wrong
    surface  the zero level extracted by marching tetrahedra
             (`tools/bench/isosurface.py`) against the exact surface:
             two-sided Hausdorff distance in cells
    time     per million nodes

    cargo run --release -p foehncraft-mesh --example sdf_benchmark -- output/sdf 48 96
    pip install numpy matplotlib
    python3 -m tools.bench.sdf output/sdf                    the table (markdown)
    python3 -m tools.bench.sdf output/sdf --json OUT         the numbers as JSON
    python3 -m tools.bench.sdf output/sdf --figure PNG       the gallery image
"""
import argparse
import json
import os
import sys

import numpy as np

from tools.bench.isosurface import hausdorff, marching_tets, rotation, sample_box, sample_sphere, sdf_box

BOX_HALF = np.array([0.55, 0.38, 0.27])
BOX_ROT = rotation(0.4, 0.3)
TORUS = (0.9, 0.35)


def sdf_torus(x):
    q = np.linalg.norm(x[..., :2], axis=-1) - TORUS[0]
    return np.sqrt(q * q + x[..., 2] ** 2) - TORUS[1]


def sample_torus(n):
    k = np.arange(n) + 0.5
    u = 2 * np.pi * k / n
    v = 2 * np.pi * ((k * 0.6180339887) % 1.0)
    big, small = TORUS
    return np.stack([(big + small * np.cos(v)) * np.cos(u), (big + small * np.cos(v)) * np.sin(u),
                     small * np.sin(v)], -1)


def exact(shape):
    if shape == "torus":
        return sdf_torus, sample_torus
    if shape == "box":
        return (lambda x: sdf_box(x, BOX_HALF, BOX_ROT)), (lambda n: sample_box(BOX_HALF, BOX_ROT, n))
    return (lambda x: np.linalg.norm(x, axis=-1) - 1.0), (lambda n: sample_sphere(1.0, n))


def load(run, shape, n, method):
    m = n + 1
    v = np.fromfile(os.path.join(run, f"{shape}_{n}_{method}.f64"), dtype="<f8")
    return v.reshape(m, m, m).transpose(2, 1, 0)  # [i, j, k] = (x, y, z)


def evaluate(run, extract_max=48):
    bench = json.load(open(os.path.join(run, "bench.json")))
    lo, hi = bench["lower"], bench["upper"]
    rows = []
    for r in bench["runs"]:
        shape, n = r["shape"], r["n"]
        h = (hi - lo) / n
        x = np.linspace(lo, hi, n + 1)
        pts = np.stack(np.meshgrid(x, x, x, indexing="ij"), -1)
        d, sample = exact(shape)
        ref = d(pts)
        near = np.abs(ref) < 0.1
        far = np.abs(ref) > h
        for method in ("ours", "sweep", "heat"):
            f = load(run, shape, n, method)
            err = np.abs(f - ref)
            row = dict(
                shape=shape, n=n, method=method,
                near_max=float(err[near].max() / h), near_rms=float(np.sqrt((err[near] ** 2).mean()) / h),
                all_max=float(err.max() / h),
                sign_errors=int((np.sign(f[far]) != np.sign(ref[far])).sum()),
                sign_fraction=float((np.sign(f[far]) != np.sign(ref[far])).mean()),
                seconds=r[f"{method}_s"], s_per_mnode=r[f"{method}_s"] / r["nodes"] * 1e6,
            )
            if n <= extract_max:
                V, T = marching_tets(pts, f)
                row["triangles"] = int(len(T))
                row["hausdorff"] = float(hausdorff(V, T, d, sample) / h) if len(T) else None
            rows.append(row)
    return bench, rows


def table(bench, rows):
    out = [f"Grid over [{bench['lower']}, {bench['upper']}]³, {bench['threads']} threads; errors in cells.", "",
           "| Shape | n | Method | Near max | Near RMS | All max | Wrong signs | Zero level Hausdorff | s per M nodes |",
           "|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        hd = r.get("hausdorff")
        out.append(
            f"| {r['shape']} | {r['n']} | {r['method']} | {r['near_max']:.3f} | {r['near_rms']:.4f} | "
            f"{r['all_max']:.2f} | {r['sign_errors']} ({100 * r['sign_fraction']:.2f} %) | "
            f"{'' if hd is None else f'{hd:.3f}'} | {r['s_per_mnode']:.3f} |")
    return "\n".join(out)


def figure(run, rows, path, n):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from tools.plot.colormaps import UNIFORM

    bench = json.load(open(os.path.join(run, "bench.json")))
    lo, hi = bench["lower"], bench["upper"]
    shapes = ["box", "holed", "soup"]
    methods = ["ours", "sweep", "heat"]
    fig, axes = plt.subplots(len(shapes), len(methods), figsize=(12, 11.4), constrained_layout=True)
    x = np.linspace(lo, hi, n + 1)
    for a, shape in enumerate(shapes):
        d, _ = exact(shape)
        X, Z = np.meshgrid(x, x, indexing="ij")
        ref = d(np.stack([X, np.zeros_like(X), Z], -1))
        for b, method in enumerate(methods):
            f = load(run, shape, n, method)[:, n // 2, :]
            ax = axes[a, b]
            im = ax.imshow(f.T, origin="lower", extent=(lo, hi, lo, hi), cmap=UNIFORM, vmin=-1.2, vmax=1.2)
            ax.contour(x, x, f.T, levels=np.arange(-1.0, 1.01, 0.2), colors="k", linewidths=0.4)
            ax.contour(x, x, ref.T, levels=[0.0], colors="w", linewidths=1.6, linestyles="--")
            ax.contour(x, x, f.T, levels=[0.0], colors="#ff2fd0", linewidths=1.4)
            r = next(r for r in rows if r["shape"] == shape and r["n"] == n and r["method"] == method)
            hd = r.get("hausdorff")
            ax.set_title(f"{shape}, {method}: near max {r['near_max']:.2f} cells,\n"
                         f"{100 * r['sign_fraction']:.1f} % wrong signs"
                         + ("" if hd is None else f", Hausdorff {hd:.2f}"), fontsize=10)
            ax.set_xticks([])
            ax.set_yticks([])
    fig.colorbar(im, ax=axes, shrink=0.6, label="signed distance φ (m); dashed white: exact surface, magenta: φ = 0")
    fig.savefig(path, dpi=100)
    print("wrote", path)


def main(argv):
    p = argparse.ArgumentParser()
    p.add_argument("run")
    p.add_argument("--json")
    p.add_argument("--figure")
    p.add_argument("--figure-n", type=int, default=96)
    a = p.parse_args(argv)
    bench, rows = evaluate(a.run)
    print(table(bench, rows))
    if a.json:
        json.dump(rows, open(a.json, "w"), indent=1)
    if a.figure:
        figure(a.run, rows, a.figure, a.figure_n)


if __name__ == "__main__":
    main(sys.argv[1:])
