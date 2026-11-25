# L17 stage A0, option B: compiler reverse mode on the momentum face kernel

Time-boxed spike of the adjoint plan
([step 27, section 2, option B](../../../docs/requirements/step27-adjoint-solver.md)):
the reverse mode of the compiler (Enzyme through rustc's
`#[autodiff_reverse]`, Moses & Churavy 2020) applied to the unchanged
kernel body `foehncraft_kernels::flow::momentum_face` on one interior face,
compared with option A (`foehncraft_kernels::reverse`, local `Dual<26>`
Jacobian, then Jᵀ ȳ). Not part of the workspace and not run in CI: it needs
a nightly toolchain, which the project allows for this spike only
(decision Q2).

```sh
rustup toolchain install nightly --profile minimal
rustup component add enzyme --toolchain nightly
cd tools/spikes/l17-enzyme-momentum-face
N=$(rustc +nightly --print sysroot)
LD_LIBRARY_PATH=$N/lib cargo +nightly rustc --release -- -Zautodiff=Enable
./target/release/l17-enzyme-momentum-face
```

`RUSTFLAGS=-Zautodiff=Enable cargo +nightly build` does not work (cargo's
target probe fails to load the Enzyme library); `cargo rustc` passes the
flag to the final crate only, which is enough. `LD_LIBRARY_PATH` is needed
because the shipped `libEnzyme-23.so` does not find rustc's LLVM library
by itself.

Result of 2026-11-24 (rustc 1.102.0-nightly 2026-11-23, LLVM 23.1.6,
one core of the cloud container): values bitwise equal, Jᵀ ȳ equal to option
A to 4.5e-16 (relative), 22 ns per primal face evaluation of the local copy,
77 ns per reverse evaluation (3.5 ×), 1 244 ns for option A recomputing
its block (56 ×). The measurement note is in
[docs/verification/adjoint-a0.md](../../../docs/verification/adjoint-a0.md).
