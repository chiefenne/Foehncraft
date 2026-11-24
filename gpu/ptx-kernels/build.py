"""Builds the single-source kernels (G1, stages 2, 3 and 5) as a PTX module
for NVIDIA GPUs.

    python3 gpu/ptx-kernels/build.py            # PTX for sm_70 (Volta), JIT-compiled for newer GPUs
    python3 gpu/ptx-kernels/build.py sm_86      # PTX for a given architecture (rebuilds `core`)

Needs rustup; the nightly toolchain, the nvptx64-nvidia-cuda target and the
llvm-bitcode-linker component are installed on first use
(rust-toolchain.toml). The rustc NVPTX target supports sm_70 and newer.
CARGO_FEATURES passes extra cargo arguments (for example
`--features adjoint`). Output:
target/nvptx64-nvidia-cuda/release/foehncraft_ptx_kernels.ptx in the
repository root (the default `run.rust_ptx` of the G1 case files).
"""

import os
import shlex
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PTX = ROOT / "target/nvptx64-nvidia-cuda/release/foehncraft_ptx_kernels.ptx"
ENTRIES = """
g1_face_flux g1_cell_gather g1_spmv_csr g1_spmv_sell
fv_gradient fv_limiter fv_faces fv_assemble fv_residual
blas_block_dot blas_block_min_max blas_update_p blas_sub_scaled blas_update_x blas_scale
blas_axpy blas_xpby
amg_smooth_first_f32 amg_smooth_step_f32 amg_add_f32 amg_residual_f32 amg_spmv_f32
amg_prolong_f32 amg_dense_f32 amg_to_f32 amg_from_f32
amg_smooth_first_f64 amg_smooth_step_f64 amg_add_f64 amg_residual_f64 amg_spmv_f64
amg_prolong_f64 amg_dense_f64 amg_to_f64 amg_from_f64
flow_lincomb flow_scaled_difference flow_fill flow_subtract flow_scale_entry
flow_block_sum_sq flow_block_min_max flow_boundary_constants flow_boundary_values
flow_gradient flow_limiter flow_mass_flux flow_momentum_faces flow_momentum_rows
flow_momentum_operator flow_residual_scale flow_set_diagonal flow_diagonal_inverse
flow_net_outflow flow_throughflow flow_courant flow_pressure_faces
flow_pressure_rows flow_correct_fluxes flow_correct_velocity
flow_neighbour_correction flow_delta_flux flow_transpose flow_matrix_product
flow_transposed_product
particle_step
""".split()


def main(argv) -> int:
    arch = argv[0] if argv else ""
    extra = shlex.split(os.environ.get("CARGO_FEATURES", ""))
    cmd = ["cargo", "build", "--release", "--target", "nvptx64-nvidia-cuda",
           "--target-dir", str(ROOT / "target")]
    env = dict(os.environ)
    if arch and arch != "sm_70":
        subprocess.run(["rustup", "component", "add", "rust-src", "--toolchain", "nightly"],
                       cwd=HERE, check=True, stdout=subprocess.DEVNULL)
        env["RUSTFLAGS"] = f"-C target-cpu={arch}"
        cmd.append("-Zbuild-std=core")
    if subprocess.run(cmd + extra, cwd=HERE, env=env).returncode != 0:
        return 1
    text = PTX.read_text(encoding="utf-8", errors="replace")
    missing = [k for k in ENTRIES if f".entry {k}" not in text]
    if missing:
        print(f"missing entry points in {PTX}: {' '.join(missing)}", file=sys.stderr)
        return 1
    target = next((line for line in text.splitlines() if line.startswith(".target")), "")
    print(f"wrote {PTX} ({target})")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
