#!/usr/bin/env python3
"""Chemistry in the flow cells (steam furnace SF2b) on a cloud GPU (Google
Colab, Kaggle, any Linux machine with an NVIDIA GPU of compute capability
7.0 or newer and its driver). Run from the repository root:

    python3 -m tools.bench.chemistry_cloud

Installs Rust if needed, builds the single-source PTX kernels for the GPU
found and the example `chemistry_gpu` of foehncraft-cuda, which integrates
a row of 4096 cells (stoichiometric methane and air from 900 to 2400 K,
fresh, half and fully burnt) over one flow time step once on the CPU
threads and once on the GPU, for GRI-Mech 3.0 and Jones & Lindstedt, and
compares the end states and step counts bitwise, with the time per cell
on both.

Writes output/chemistry_gpu/report.md (also printed) and exits with 1 when
a GPU result differs from its CPU result. Takes about 10 minutes on a free
Colab T4 (most of it compiling).
"""
import subprocess
import sys

from tools.bench.stage5_cloud import ROOT, gpu, rust, sh

OUT = ROOT / "output" / "chemistry_gpu"


def main():
    name, sm = gpu()
    rust()
    print(f"== Single-source kernels -> PTX for sm_{sm}", flush=True)
    sh([sys.executable, "gpu/ptx-kernels/build.py", f"sm_{sm}"], check=True)
    print("== Chemistry example with CUDA", flush=True)
    sh(["cargo", "build", "--release", "--locked", "-p", "foehncraft-cuda", "--features", "cuda",
        "--example", "chemistry_gpu"], check=True)
    exe = ROOT / "target/release/examples/chemistry_gpu"
    ptx = ROOT / "target/nvptx64-nvidia-cuda/release/foehncraft_ptx_kernels.ptx"
    r = subprocess.run([str(exe), str(ptx)], cwd=ROOT, capture_output=True, text=True, check=False)
    sys.stderr.write(r.stderr)
    rep = f"# Chemistry SF2b: CPU against GPU\n\nGPU: {name} (sm_{sm}).\n\n" + r.stdout
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "report.md").write_text(rep)
    print("\n" + rep)
    return r.returncode


if __name__ == "__main__":
    sys.exit(main())
