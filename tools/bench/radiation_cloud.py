#!/usr/bin/env python3
"""Radiation R1a on a cloud GPU (Google Colab, Kaggle, any Linux machine
with an NVIDIA GPU of compute capability 7.0 or newer and its driver). Run
from the repository root:

    python3 -m tools.bench.radiation_cloud

Installs Rust if needed, builds the single-source PTX kernels for the GPU
found and the example `radiation_gpu` of foehncraft-cuda, which runs every
radiation case once with the bundles on the CPU threads and once on the
GPU and compares the results bitwise (D259): the slab (case 1), the
isothermal enclosure (case 2), the parallel plates (case 3) and a larger
case with every wall kind on about 30 000 polyhedra, with the time per
bundle on both.

Writes output/radiation_gpu/report.md (also printed) and exits with 1 when
a GPU result differs from its CPU result. Takes about 10 minutes on a free
Colab T4 (most of it compiling).
"""
import subprocess
import sys

from tools.bench.stage5_cloud import ROOT, gpu, rust, sh

OUT = ROOT / "output" / "radiation_gpu"


def main():
    name, sm = gpu()
    rust()
    print(f"== Single-source kernels -> PTX for sm_{sm}", flush=True)
    sh([sys.executable, "gpu/ptx-kernels/build.py", f"sm_{sm}"], check=True)
    print("== Radiation example with CUDA", flush=True)
    sh(["cargo", "build", "--release", "--locked", "-p", "foehncraft-cuda", "--features", "cuda",
        "--example", "radiation_gpu"], check=True)
    exe = ROOT / "target/release/examples/radiation_gpu"
    ptx = ROOT / "target/nvptx64-nvidia-cuda/release/foehncraft_ptx_kernels.ptx"
    r = subprocess.run([str(exe), str(ptx)], cwd=ROOT, capture_output=True, text=True, check=False)
    sys.stderr.write(r.stderr)
    rep = f"# Radiation R1a: CPU against GPU\n\nGPU: {name} (sm_{sm}).\n\n" + r.stdout
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "report.md").write_text(rep)
    print("\n" + rep)
    return r.returncode


if __name__ == "__main__":
    sys.exit(main())
