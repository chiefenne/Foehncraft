#!/usr/bin/env python3
"""The deforming mesh (stage MS0, D91) on a cloud GPU (Google Colab,
Kaggle, any Linux machine with an NVIDIA GPU of compute capability 7.0 or
newer and its driver). Run from the repository root:

    python3 -m tools.bench.motion_cloud

Installs Rust if needed, builds the single-source PTX kernels for the GPU
found and the CLI with the `cuda` feature, then runs the deforming-mesh
cases twice, with the flow loop on the CPU (`--backend cpu`) and on the
GPU (`--backend cuda`), on one rank (ranks run on the CPU device):

- motion_free_stream: a uniform flow on a deforming mesh stays uniform
  (MS0.1); the GPU run must keep it to round-off as well;
- motion_channel_ranks: the breathing channel with a dye stripe, and its
  continuation from restart_000075.h5 on the GPU (MS0 part 2);
- motion_piston_adiabatic: the L22 piston with the variable density.

Since MPI step W2 the CPU sums exactly and the GPU keeps its block sums
(owner decision of 2026-11-23), so the two runs agree to rounding, not
bitwise: the report gives the largest difference of every field of the
last field file relative to the field's range, and the physics checks of
each case on the GPU. Exits with 1 when a run fails, a check fails or a
difference exceeds 1e-9. Writes output/motion_gpu/report.md (also
printed). Takes about 15 minutes on a free Colab T4 (most of it the build).
"""
import json
import pathlib
import shutil
import sys

import h5py
import numpy as np

from tools.bench.stage5_cloud import ROOT, gpu, rust, sh

OUT = ROOT / "output" / "motion_gpu"
F = "examples/cases/flow/"
CASES = ["motion_free_stream", "motion_channel_ranks", "motion_piston_adiabatic"]
LIMIT = 1e-9


def fields(d):
    """The cell fields of the last field file of run directory d."""
    out = {}
    with h5py.File(d / "fields.vtkhdf", "r") as f:
        g = f["VTKHDF"]["internal"]["CellData"]
        for k in g:
            if k in ("GlobalCellIds", "OriginalCellIds", "vtkGhostType", "zone"):
                continue
            out[k] = np.asarray(g[k][()], dtype=float)
    return out


def worst(a, b):
    """The largest difference of a field relative to its range, and its name."""
    w, name = 0.0, ""
    for k, va in a.items():
        vb = b.get(k)
        if vb is None or va.shape != vb.shape:
            return float("inf"), k
        scale = max(np.ptp(va), np.abs(va).max(), 1e-300)
        d = float(np.abs(va - vb).max()) / scale
        if d > w:
            w, name = d, k
    return w, name


def solve(d):
    return json.loads((d / "result.json").read_text())["solves"][-1]


def check(case, d):
    """The physics check of a case on the GPU: (ok, text)."""
    s = solve(d)
    m = s["transient"]["mesh_motion"]
    text = f"closure {m['closure_max']:.1e}"
    ok = m["closure_max"] < 1e-12
    if case == "motion_free_stream":
        u = fields(d)
        vel = [k for k in u if k.startswith("velocity_")]
        spread = max(float(np.ptp(u[k])) for k in vel) if vel else float("inf")
        text += f", velocity spread {spread:.1e} m/s"
        ok &= spread < 1e-9
    if case == "motion_piston_adiabatic":
        r, cp, p_op, v0 = 8.314462618 / 0.0289647, 1005.0, 101325.0, 0.02 * 0.04 * 0.002
        gamma = cp / (cp - r)
        exact = p_op * (v0 / m["volume"]) ** gamma
        p0 = s["density"]["p0"]
        text += f", p0 {p0:.6e} Pa against {exact:.6e} (p0 V^gamma)"
        ok &= abs(p0 / exact - 1) < 1e-3
    return ok, text


def main():
    name, sm = gpu()
    rust()
    print(f"== Single-source kernels -> PTX for sm_{sm}", flush=True)
    sh([sys.executable, "gpu/ptx-kernels/build.py", f"sm_{sm}"], check=True)
    print("== CLI with CUDA", flush=True)
    sh(["cargo", "build", "--release", "--locked", "-p", "foehncraft-cli", "--features", "cuda"], check=True)
    foehn = str(ROOT / "target/release/foehn")
    shutil.rmtree(OUT, ignore_errors=True)
    rows, status = [], 0
    runs = [(c, None) for c in CASES] + [("motion_channel_ranks", "restart")]
    for case, mode in runs:
        tag = case + ("_restart" if mode else "")
        ok = True
        for be in ("cpu", "cuda"):
            d = OUT / tag / be
            cmd = [foehn, "flow", F + case + ".toml", "--backend", be, "--ranks", "1", "--out", str(d)]
            if mode:
                cmd += ["--restart", str(OUT / case / be / "restart_000075.h5")]
            if sh(cmd).returncode != 0:
                print(f"!! {tag} on {be} failed", flush=True)
                ok = False
        if not ok:
            rows.append(f"| {tag} | failed | | |")
            status = 1
            continue
        w, field = worst(fields(OUT / tag / "cpu"), fields(OUT / tag / "cuda"))
        good, text = check(case, OUT / tag / "cuda")
        if mode:
            # the continued GPU run ends where the GPU run without a break ends
            c, _ = worst(fields(OUT / case / "cuda"), fields(OUT / tag / "cuda"))
            text += f"; continued against unbroken GPU run {c:.1e}"
            good &= c == 0.0
        good &= w <= LIMIT
        status |= 0 if good else 1
        rows.append(f"| {tag} | {w:.1e} ({field}) | {text} | {'yes' if good else '**no**'} |")
    rep = (
        f"# Deforming mesh (MS0) on the GPU\n\nGPU: {name} (sm_{sm}).\n\n"
        "| case | largest CPU-GPU difference / range | GPU checks | passed |\n|---|---|---|---|\n"
        + "\n".join(rows) + "\n"
    )
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "report.md").write_text(rep)
    print("\n" + rep)
    print("all deforming-mesh GPU runs pass" if status == 0 else "FAILED: see above")
    return status


if __name__ == "__main__":
    sys.exit(main())
