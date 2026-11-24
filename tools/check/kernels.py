"""R10: the kernels stay reproducible (D31) and on the float64 path (D7).

    python -m tools.check.kernels [PTX]   # default: the module of python3 gpu/ptx-kernels/build.py

Fails when
- the kernel sources use atomics (no floating-point atomics, no atomic
  scatter: every output is gathered by one index), except the integer
  tallies of the Monte Carlo radiation (D269: `gpu/ptx-kernels/src/radiation.rs`);
- the kernel sources use f32 outside the multigrid modules (`amg.rs`: the
  float32 hierarchy of D37; everything else stays on the float64 path);
- the PTX module contains atomic or reduction instructions (atom., red.)
  other than 64-bit integer adds in the entries of INTEGER_ATOMICS (D269:
  a sum of integers is the same in any order, so it stays reproducible),
  fused multiply-adds (fma.rn.f64, fma.rn.f32: results would differ from
  the CPU), approximate division or reciprocals, or flush-to-zero f32
  arithmetic (.ftz).
"""

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SOURCES = [ROOT / "crates/foehncraft-kernels/src", ROOT / "gpu/ptx-kernels/src"]
DEFAULT_PTX = ROOT / "target/nvptx64-nvidia-cuda/release/foehncraft_ptx_kernels.ptx"
ATOMIC = re.compile(r"Atomic|atomic_|core::intrinsics::atomic")
F32 = re.compile(r"\bf32\b")
# comment lines and the Real impl of f32 are not uses
F32_OK = re.compile(r"^\s*//|impl Real for f32|//!|///")
PTX_ATOMIC = re.compile(r"^\s*(atom|red)\.")
# D269: entries that may add 64-bit integers atomically, and the files
# that may write the atomics
INTEGER_ATOMICS = {"radiation_transport"}
ATOMIC_SOURCES = {"gpu/ptx-kernels/src/radiation.rs"}
PTX_INTEGER_ADD = re.compile(r"^\s*(atom|red)(\.relaxed)?(\.(gpu|sys|cta))?(\.global)?\.add\.[us]64\s")
PTX_FUNCTION = re.compile(r"^\.(visible\s+)?(\.?entry|\.?func)\b.*?([A-Za-z_][A-Za-z0-9_$]*)\s*\(")
PTX_INEXACT = re.compile(r"fma\.rn\.f(64|32)|div\.(approx|full)\.f(64|32)|rcp\.approx|\.ftz\.")


def grep(pattern, files, skip=None):
    hits = []
    for f in files:
        for n, line in enumerate(f.read_text(encoding="utf-8", errors="replace").splitlines(), 1):
            if pattern.search(line) and not (skip and skip(f, line)):
                hits.append(f"{f.relative_to(ROOT) if f.is_relative_to(ROOT) else f}:{n}:{line}")
    return hits


def main(argv) -> int:
    if argv and argv[0] in ("-h", "--help"):
        print(__doc__)
        return 0
    ptx = Path(argv[0]) if argv else DEFAULT_PTX
    src = sorted(p for d in SOURCES for p in d.rglob("*") if p.is_file())
    fail = False
    hits = grep(ATOMIC, src, lambda f, line: f.relative_to(ROOT).as_posix() in ATOMIC_SOURCES)
    if hits:
        print("\n".join(hits))
        print("R10: atomics in the kernel sources (above)", file=sys.stderr)
        fail = True
    hits = grep(F32, src, lambda f, line: f.name == "amg.rs" or F32_OK.search(line))
    if hits:
        print("\n".join(hits))
        print("R10: f32 outside the multigrid kernels (above)", file=sys.stderr)
        fail = True
    if ptx.is_file():
        hits = []
        function = None
        for n, line in enumerate(ptx.read_text(encoding="utf-8", errors="replace").splitlines(), 1):
            m = PTX_FUNCTION.match(line)
            if m:
                function = m.group(3)
            if PTX_ATOMIC.search(line) and not (
                function in INTEGER_ATOMICS and PTX_INTEGER_ADD.search(line)
            ):
                hits.append(f"{ptx}:{n}:{line}")
        if hits:
            print("\n".join(hits))
            print(f"R10: atomic instructions in {ptx} (above)", file=sys.stderr)
            fail = True
        hits = grep(PTX_INEXACT, [ptx])
        if hits:
            print("\n".join(hits))
            print(f"R10: contracted, approximate or flush-to-zero arithmetic in {ptx} (above)", file=sys.stderr)
            fail = True
        print(f"PTX checked: {ptx}")
    else:
        print(f"no PTX module at {ptx} (build it with python3 gpu/ptx-kernels/build.py); sources checked only")
    if not fail:
        print("R10 kernel checks passed")
    return 1 if fail else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
