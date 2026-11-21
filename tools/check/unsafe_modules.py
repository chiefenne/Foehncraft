#!/usr/bin/env python3
"""D202: `unsafe` only in the named GPU and SIMD modules.

    python3 -m tools.check.unsafe_modules

The workspace denies `unsafe_code` (Cargo.toml, `[workspace.lints.rust]`)
and clippy's `undocumented_unsafe_blocks`, so the compiler already rejects
an unsafe block outside an allowed module, and one without a `// SAFETY:`
comment. A module is allowed by `#![allow(unsafe_code)]` at its top; this
check keeps that list short and named: it fails (exit code 1) when

- a Rust file outside ALLOWED allows `unsafe_code` (inner or outer
  attribute, `allow` or `expect`), or
- a file of ALLOWED no longer allows it (the entry is stale), or
- a Rust file of the workspace crates outside ALLOWED (and outside a crate
  allowed as a whole) contains the keyword `unsafe` in code.

A new unsafe module needs the maintainer's go and one line here, with the
reason.
"""
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]

# file (or crate directory, ending in "/") -> why it needs unsafe code
ALLOWED = {
    # GPU: the kernel bodies are one source for the CPU and the NVIDIA GPU;
    # their loads are unchecked behind `unsafe fn` (crate docs, "Safety")
    "crates/foehncraft-kernels/": "kernel bodies, single source for CPU and GPU",
    # GPU: driver calls, device memory, kernel launches (D60)
    "crates/foehncraft-cuda/": "CUDA backend",
    # SIMD and CPU loops of the kernel bodies (CSR and SELL-C rows)
    "crates/foehncraft-fv/src/backend/cpu.rs": "CPU backend loops over the kernel bodies",
    "crates/foehncraft-bench/src/backend/cpu.rs": "CPU backend of gate G1",
    "crates/foehncraft-bench/src/workload.rs": "G1 reference flux from the face kernel",
    "crates/foehncraft-bench/src/matrix.rs": "tests of the CSR and SELL-C row kernels",
    "crates/foehncraft-fv/src/amg/csr.rs": "test of the SELL-C row kernel",
    "crates/foehncraft-flow/src/device/cpu.rs": "host executor: device arrays, launches",
    # foreign function interface
    "crates/foehncraft-fv/src/hypre.rs": "C interface of the hypre shim",
    # tests that call the kernel bodies one face or cell at a time
    "crates/foehncraft-flow/tests/adjoint.rs": "kernel bodies under test",
    "crates/foehncraft-flow/tests/quality.rs": "kernel bodies under test",
    "crates/foehncraft-flow/tests/symmetry_stress.rs": "kernel bodies under test",
}
SKIP_DIRS = {"target", ".git", "node_modules", "__pycache__"}
ALLOW = re.compile(r"#!?\[\s*(allow|expect)\s*\([^)]*\bunsafe_code\b")
KEYWORD = re.compile(r"\bunsafe\b")


def allowed(rel):
    return any(rel == a or (a.endswith("/") and rel.startswith(a)) for a in ALLOWED)


def code(line):
    """The code of a line: comments and string contents are not code."""
    s = line.split("//", 1)[0]
    return re.sub(r'"([^"\\]|\\.)*"', '""', s)


def main():
    errors = []
    seen = set()
    for path in sorted((ROOT / "crates").rglob("*.rs")):
        if SKIP_DIRS & set(path.parts):
            continue
        rel = path.relative_to(ROOT).as_posix()
        text = path.read_text(encoding="utf-8")
        has_allow = any(ALLOW.search(code(l)) for l in text.splitlines())
        if rel in ALLOWED:
            seen.add(rel)
            if not has_allow:
                errors.append(f"{rel}: listed in ALLOWED but has no #![allow(unsafe_code)]")
            continue
        if allowed(rel):
            continue
        if has_allow:
            errors.append(f"{rel}: allows unsafe_code but is not a named unsafe module")
        for no, line in enumerate(text.splitlines(), 1):
            if KEYWORD.search(code(line)):
                errors.append(f"{rel}:{no}: `unsafe` outside the named modules")
    for a in ALLOWED:
        if not a.endswith("/") and a not in seen:
            errors.append(f"{a}: listed in ALLOWED but does not exist")
    for e in errors:
        print(e)
    if errors:
        print(f"{len(errors)} finding(s); see tools/check/unsafe_modules.py (D202)")
        return 1
    print(f"unsafe code only in the {len(ALLOWED)} named modules (D202)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
