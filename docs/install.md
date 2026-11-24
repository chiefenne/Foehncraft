# Install and update

Foehncraft is built from source, also release 1.0.1 (tag `v1.0.1`);
release binaries for Linux, macOS and Windows with `SHA256SUMS.txt` and
`foehn doctor` are planned (D203) but not part of a release yet. See also
[Connect an assistant](agents.md) and [Troubleshooting](troubleshooting.md).

## What you need

- Rust 1.89 or newer (`rustup`, <https://rustup.rs>).
- A C compiler and CMake: HDF5 (VTKHDF output, restart files) and METIS
  (partitioning) are built from source on the first build.
- Git.
- Optional: an NVIDIA driver for the GPU backend (`--features cuda`; no CUDA
  toolkit is needed to build), Python 3.10 or newer for the Python package.

## Install

```
git clone https://github.com/chiefenne/Foehncraft.git
cd Foehncraft
git checkout v1.0.1          # optional: the release instead of main
cargo build --release -p foehncraft-cli
target/release/foehn info
```

The first build takes several minutes (HDF5 and METIS); later builds are
fast. To call `foehn` from anywhere, copy the binary onto your `PATH`:

```
# Linux and macOS
install -m 755 target/release/foehn ~/.local/bin/foehn
# Windows (PowerShell)
Copy-Item target\release\foehn.exe $env:USERPROFILE\bin\
```

or install it with Cargo: `cargo install --path crates/foehncraft-cli`.

With the GPU backend: `cargo build --release -p foehncraft-cli --features cuda`.
Without the C libraries (no VTKHDF, no METIS):
`cargo build --release -p foehncraft-cli --no-default-features`.

The Python package: `pip install ./python` ([recipe 33](recipes/python.md#python-api)).

## Check the installation

```
foehn info                                    # version, commit, features, platform
foehn commands | head                         # the tool catalogue
foehn flow examples/cases/flow/cavity_re100.toml --out output/cavity
foehn mcp --read-root . --write-root output < /dev/null   # starts and ends at once
```

## Update

```
cd Foehncraft
git checkout main && git pull
cargo build --release -p foehncraft-cli
foehn info                                    # the new commit
```

Then, if you copied the binary, copy it again. An assistant host keeps
the old server running until you restart the host or the server.

After an update:

- `foehn case upgrade FILE` brings case files of an older case format to
  the current one (`--write` rewrites them, keeping their comments).
- A journal of an older version replays with `foehn run --script`; a call
  that ends otherwise is reported (`TOOL005`).
- The [changelog](../CHANGELOG.md) lists what changed; the "Unreleased"
  part is on `main` but not yet in a release.
