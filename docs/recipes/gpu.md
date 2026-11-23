# GPU runs

Part of the [Recipes](../recipes.md): scalar operators, pressure solver and flow cases on an NVIDIA GPU (also a free Colab T4).

| # | Recipe |
|---|---|
| 14 | [Run the scalar operators on an NVIDIA GPU](#scalar-gpu) (also: free Colab T4) |
| 16 | [Run the pressure solver on an NVIDIA GPU against hypre](#poisson-gpu) (also: free Colab T4) |
| 24 | [Run a flow case on an NVIDIA GPU and compare it with the CPU](#flow-gpu) (also: free Colab T4) |

---

<a id="scalar-gpu"></a>
## 14. Run the scalar operators on an NVIDIA GPU

Without a GPU of your own: open
[`tools/notebooks/stage2_colab.ipynb`](https://colab.research.google.com/github/chiefenne/Foehncraft/blob/main/tools/notebooks/stage2_colab.ipynb)
in Google Colab, choose a T4 runtime, set `BRANCH` if the stage-2 branch is
not merged yet, and run all cells (about 15 minutes). On any Linux machine
with an NVIDIA GPU of compute capability 7.0 or newer:

```
python3 -m tools.bench.stage2_cloud
```

By hand:

```
python3 gpu/ptx-kernels/build.py sm_75                          # your GPU's sm_XY (nvidia-smi: compute_cap)
cargo build --release -p foehncraft-cli --features cuda
target/release/foehn scalar examples/cases/scalar/mms_poly_gpu.toml
target/release/foehn scalar examples/cases/scalar/oblique_step_2d.toml --backends cpu,cuda
target/release/foehn scalar examples/cases/scalar/mms_poly_box.toml --backends cuda
```

`mms_poly_gpu.toml` solves on the CPU and twice on the GPU (`run.backends =
["cpu", "cuda", "cuda"]`); the column "vs first backend" of the report says
"bitwise equal" when a GPU run reproduces the CPU run bit for bit, and the
run fails if any difference exceeds 1e-10. Since MPI step W2 (#219) the CPU sums exactly and the GPU keeps block sums, so GPU results agree with the CPU to rounding, not bitwise, until the GPU sums exactly too (v1.1, owner decision of 2026-11-23). Keys: `run.backends`,
`run.rust_ptx` (or `--rust-ptx`).

[↑ back to top](#gpu-runs)

---

<a id="poisson-gpu"></a>
## 16. Run the pressure solver on an NVIDIA GPU against hypre

Without a GPU of your own: open
[`tools/notebooks/stage3_colab.ipynb`](https://colab.research.google.com/github/chiefenne/Foehncraft/blob/main/tools/notebooks/stage3_colab.ipynb)
in Google Colab, choose a T4 runtime, set `BRANCH` if the stage-3 branch is
not merged yet, and run all cells (about 30 minutes, most of it compiling
hypre with CUDA). On any Linux machine with an NVIDIA GPU of compute
capability 7.0 or newer and the CUDA toolkit:

```
python3 -m tools.bench.stage3_cloud
```

By hand (sm_75 for a T4; `nvidia-smi --query-gpu=compute_cap --format=csv`):

```
python3 -m tools.bench.build_hypre $PWD/target/hypre-cuda75 75
python3 gpu/ptx-kernels/build.py sm_75
HYPRE_DIR=$PWD/target/hypre-cuda75 cargo build --release -p foehncraft-cli --features cuda,hypre
target/release/foehn poisson examples/cases/poisson/poisson_poly_gpu.toml
target/release/foehn poisson examples/cases/poisson/poisson_poly_box.toml --backends cuda --solvers amg,amg-double,hypre
```

`poisson_poly_gpu.toml` solves 912 673 polyhedra with the in-house solvers
on the CPU and on the GPU and with hypre on the GPU; the report gives the
setup and solve times and ns per cell and iteration of each, and the
column "Bitwise vs first backend" says whether the GPU solve equals the CPU
solve bit for bit. Keys: `study.backends`, `study.solvers`,
`study.rust_ptx` (or `--rust-ptx`), `study.repeat`.

[↑ back to top](#gpu-runs)

---

<a id="flow-gpu"></a>
## 24. Run a flow case on an NVIDIA GPU and compare it with the CPU

Every flow case runs on the GPU as it is: the whole flow loop (gradients,
limiters, momentum and pressure-correction assembly, Rhie–Chow fluxes,
corrections, the BiCGStab and CG + multigrid solves) stays on the device.
Up to v0.7 the results were bitwise those of the CPU (D31); Since MPI step W2 (#219) the CPU sums exactly and the GPU keeps block sums, so GPU results agree with the CPU to rounding, not bitwise, until the GPU sums exactly too (v1.1, owner decision of 2026-11-23). Needs an NVIDIA GPU of
compute capability 7.0 or newer with its driver; the build needs no CUDA
toolkit, the kernels need the nightly toolchain (installed by
`rust-toolchain.toml` on first use). From the repository root:

```
python3 gpu/ptx-kernels/build.py sm_75                       # PTX for your GPU (sm_75: T4, sm_86: RTX 30xx, ...)
cargo build --release -p foehncraft-cli --features cuda
target/release/foehn info                             # features: cuda ...
target/release/foehn flow examples/cases/flow/cavity_re100.toml --backend cuda
```

To check the GPU against the CPU, run the case twice and compare the output
directories (`pip install h5py numpy` once):

```
target/release/foehn flow examples/cases/flow/cylinder_re100_gpu.toml --backend cpu  --out output/gpu/cpu
target/release/foehn flow examples/cases/flow/cylinder_re100_gpu.toml --backend cuda --out output/gpu/cuda
python3 -m tools.check.compare_runs output/gpu/cpu output/gpu/cuda    # until v1.1: differences in the last bits
```

Until the GPU sums exactly (v1.1), `compare_runs` reports the CPU and GPU
outputs as different in the last bits; two CPU runs with different thread
counts stay bitwise equal.

`python3 -m tools.bench.stage5_cloud` does this for eight cases, from Couette
flow to the cavity on one million polyhedra, and writes
`output/stage5_gpu/report.md` with the times per cell and iteration.
**Without a GPU of your own:** open
[`tools/notebooks/stage5_colab.ipynb`](https://colab.research.google.com/github/chiefenne/Foehncraft/blob/main/tools/notebooks/stage5_colab.ipynb)
in Google Colab, choose a T4 runtime, set `BRANCH` if the stage-5 branch
is not merged yet, and run all cells (about 15 minutes).

**The turbulence models (stage 6, specification 6.6):**
`python3 -m tools.bench.stage6_cloud` runs every model on the TMR flat
plate and the channel, the heated channel and pipe and the coupled step,
300 outer iterations each, on the CPU and on the GPU, and compares each
pair bitwise (`output/stage6_gpu/report.md`); `--quick` runs the five
channels only. Without a GPU, `--cpu-only` compares one thread with four
on the CPU. In Colab:
[`tools/notebooks/stage6_colab.ipynb`](https://colab.research.google.com/github/chiefenne/Foehncraft/blob/main/tools/notebooks/stage6_colab.ipynb)
(T4 runtime, run all cells, about 20 minutes).

Keys: `run.backend` (`cpu`, `cuda`), `run.ptx` (the PTX module); the
command-line options `--backend` and `--ptx` override them. A build
without the `cuda` feature, or a machine without the NVIDIA driver, stops
with a message instead of running on the CPU.

[↑ back to top](#gpu-runs)
