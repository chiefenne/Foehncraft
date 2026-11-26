# Install and update

Foehncraft is built from source, also release 1.1.0 (tag `v1.1.0`). Release binaries for Linux, macOS and
Windows with `SHA256SUMS.txt` and `foehn doctor` are planned (D203) but not
part of a release yet. This page lists every tool to install on macOS,
Linux and Windows, with the commands, so that a fresh machine gets to a
working `foehn` in one pass. See also [Connect an assistant](agents.md) and
[Troubleshooting](troubleshooting.md).

Contents: [What you need](#what-you-need) ·
[macOS](#macos) · [Linux](#linux) · [Windows](#windows) ·
[Build](#build) · [Python, Gmsh and the example scripts](#python-gmsh-and-the-example-scripts) ·
[Optional: MPI and GPU](#optional-mpi-and-gpu) ·
[Check the installation](#check-the-installation) · [Update](#update)

## What you need

| Tool | Needed for | Required |
|---|---|---|
| Rust 1.89 or newer (rustup) | building Foehncraft | yes |
| C and C++ compiler | HDF5, METIS and BLAKE3, built from source on the first build | yes |
| CMake 3.10 or newer | the HDF5 and METIS builds | yes |
| Git | cloning the repository | yes |
| Python 3.11 or newer (3.12 in CI) | the Python package, the example mesh scripts, plots and studies | for scripts |
| Gmsh Python module (`pip install gmsh`) | the mesh scripts in `examples/meshes/`, `examples/applications/*/mesh.py` and some studies | for those scripts |
| Open MPI and libclang | runs on several MPI ranks (`--features mpi`) | optional |
| NVIDIA driver | the GPU backend (`--features cuda`; no CUDA toolkit needed) | optional |

You do **not** need a system HDF5, METIS or the Gmsh program: Foehncraft
builds HDF5 and METIS itself (that is what CMake is for) and reads Gmsh
`.msh` files with its own reader. Gmsh is only needed to *make* new meshes
with the scripts in `examples/`.

CI builds and tests on `ubuntu-latest`, `macos-latest` and
`windows-latest` with exactly these tools (the GitHub runners come with
CMake and the C compilers; CI adds only Open MPI for the MPI job).

## macOS

1. The Xcode Command Line Tools (C and C++ compiler, Git, libclang):

   ```
   xcode-select --install
   ```

2. [Homebrew](https://brew.sh), if you do not have it yet:

   ```
   /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
   ```

3. CMake and Python:

   ```
   brew install cmake python@3.12
   ```

4. Rust with rustup (accept the default installation), then open a new
   terminal:

   ```
   curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
   rustc --version          # 1.89 or newer; else: rustup update stable
   cmake --version
   ```

   With a Homebrew Rust (`brew install rust`) instead of rustup, check the
   version the same way.

Continue with [Build](#build).

## Linux

### Ubuntu and Debian

```
sudo apt update
sudo apt install build-essential cmake git curl python3 python3-venv python3-pip
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source "$HOME/.cargo/env"
rustc --version          # 1.89 or newer; else: rustup update stable
```

### Fedora, RHEL, Rocky and Alma

```
sudo dnf install gcc gcc-c++ make cmake git curl python3 python3-pip
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source "$HOME/.cargo/env"
rustc --version
```

(The Ubuntu commands are the ones checked on Ubuntu 24.04; the Fedora
package names are the usual equivalents.)

Continue with [Build](#build).

## Windows

Use PowerShell. The commands use `winget`, which comes with Windows 10 and
11.

1. The Visual Studio C++ Build Tools (the MSVC compiler and linker that Rust
   and the HDF5 build need):

   ```
   winget install --id Microsoft.VisualStudio.2022.BuildTools -e --override "--quiet --wait --add Microsoft.VisualStudio.Workload.VCTools --includeRecommended"
   ```

   (`rustup-init.exe` in the next step offers to install them too, if you
   skip this.)

2. Rust, CMake, Git and Python:

   ```
   winget install --id Rustlang.Rustup -e
   winget install --id Kitware.CMake -e
   winget install --id Git.Git -e
   winget install --id Python.Python.3.12 -e
   ```

3. Close PowerShell and open a new one (so that the new `PATH` applies),
   then check:

   ```
   rustc --version          # 1.89 or newer; else: rustup update stable
   cmake --version
   git --version
   py --version
   ```

With [Chocolatey](https://chocolatey.org) instead of winget:
`choco install rustup.install git python visualstudio2022buildtools visualstudio2022-workload-vctools` and
`choco install cmake --installargs 'ADD_CMAKE_TO_PATH=System'` (without
that argument CMake is not put on `PATH`).

Continue with [Build](#build). In the commands below, the program is
`target\release\foehn.exe`; PowerShell accepts the `/` paths as written.

## Build

On all three systems, in a terminal (PowerShell on Windows):

```
git clone https://github.com/chiefenne/Foehncraft.git
cd Foehncraft
git checkout v1.1.0          # optional: the release instead of main
cargo build --release --locked -p foehncraft-cli
target/release/foehn info
```

To build a release instead of `main`, check out its tag before the build
(`git tag` lists them, for example `git checkout v1.0.1`).

The first build takes several minutes (HDF5 and METIS are compiled from
source); later builds are fast. To call `foehn` from any folder, install it
into `~/.cargo/bin`, which rustup puts on your `PATH`:

```
cargo install --locked --path crates/foehncraft-cli
foehn info
```

Without the C libraries (no VTKHDF output, no restart files, no METIS
partitioning; no CMake needed):
`cargo build --release --locked -p foehncraft-cli --no-default-features`.

## Python, Gmsh and the example scripts

The scripts in `examples/` and `tools/` are Python. Make a virtual
environment once, in the `Foehncraft` folder:

```
# macOS and Linux
python3 -m venv .venv
source .venv/bin/activate

# Windows (PowerShell)
py -3.12 -m venv .venv
.venv\Scripts\Activate.ps1
```

Then install what the scripts use:

```
pip install gmsh numpy matplotlib h5py      # mesh scripts, plots, fields
pip install ./python                        # the foehncraft Python package (builds the Rust module)
```

- **Gmsh** comes as the pip package `gmsh` (the Gmsh library with its Python
  API); no separate Gmsh program is needed. Check it with
  `python -c "import gmsh; gmsh.initialize(); print(gmsh.__version__); gmsh.finalize()"`.
- **On Linux** the `gmsh` package needs some OpenGL and X11 libraries even
  without a desktop (otherwise `import gmsh` stops with
  `OSError: libGLU.so.1: cannot open shared object file`):

  ```
  # Ubuntu and Debian
  sudo apt install libglu1-mesa libxcursor1 libxft2 libxinerama1 libgomp1
  # Fedora
  sudo dnf install mesa-libGLU libXcursor libXft libXinerama libgomp
  ```

- `pip install ./python` builds the native module with maturin, so it needs
  the Rust and CMake setup above; extras: `pip install "./python[fields]"`
  (h5py, numpy), `"./python[study]"` (scipy). See
  [recipe 33](recipes/python.md#python-api).
- Some scripts name further packages in their first lines (for example
  `scipy`, `pyvista`, `pillow`, `sympy`); install them the same way.

A first mesh, to check Gmsh:

```
python examples/meshes/mixer.py --out output/mixer_mesh
```

## Optional: MPI and GPU

**MPI** (`foehn run --ranks N`, several processes) needs Open MPI and
libclang (for the bindings):

```
# macOS (libclang comes with the Xcode Command Line Tools)
brew install open-mpi
# Ubuntu and Debian
sudo apt install libopenmpi-dev openmpi-bin libclang-dev
# Fedora (then load the module: module load mpi/openmpi-x86_64)
sudo dnf install openmpi-devel clang-devel
```

```
cargo build --release --locked -p foehncraft-cli --features mpi
mpirun -n 4 target/release/foehn flow examples/cases/flow/ranks_run_channel.toml --out output/ranks
```

Add `--oversubscribe` to `mpirun` when the machine has fewer
cores than ranks. CI tests MPI on Linux and macOS; on Windows MPI is not built or tested.

**GPU**: `cargo build --release --locked -p foehncraft-cli --features cuda`
needs an NVIDIA driver at run time and no CUDA toolkit to build. The GPU
kernels themselves need nightly Rust (see `gpu/ptx-kernels`); the Colab
notebooks in `tools/` show a working setup.

## Check the installation

```
foehn info                                    # version, commit, features, platform
foehn commands | head                         # the tool catalogue
foehn flow examples/cases/flow/cavity_re100.toml --out output/cavity
foehn mcp --read-root . --write-root output < /dev/null   # starts and ends at once
```

(Without `cargo install`, write `target/release/foehn` for `foehn`. On
Windows, leave out `| head` and `< /dev/null`.)

## Update

```
git checkout main && git pull --ff-only
cargo build --release --locked -p foehncraft-cli
foehn info                                    # the new commit
```

If you installed it with `cargo install`, run that again, and
`pip install ./python` again for the Python package. An assistant host
keeps the old server running until you restart the host or the server.
After `rustup update stable` nothing else is needed.

After an update:

- `foehn case upgrade FILE` brings case files of an older case format to
  the current one (`--write` rewrites them, keeping their comments).
- A journal of an older version replays with `foehn run --script`; a call
  that ends otherwise is reported (`TOOL005`).
- The [changelog](../CHANGELOG.md) lists what changed; the "Unreleased"
  part is on `main` but not yet in a release.

## Problems on a first build

| Message | Cause and fix |
|---|---|
| `failed to execute command: ... "cmake"` or `is cmake not installed?` | CMake is missing or not on `PATH`: install it as above, open a new terminal. |
| `linker 'cc' not found` (Linux), `link.exe not found` (Windows) | No C compiler: `build-essential` / `gcc`, or the Visual Studio C++ Build Tools. |
| `error: package ... requires rustc 1.89` | `rustup update stable`. |
| `Unable to find libclang` | Only with `--features mpi`: install `libclang-dev` (Ubuntu) or `clang-devel` (Fedora). |
| `OSError: libGLU.so.1` on `import gmsh` | The Linux libraries of the Gmsh section above. |
