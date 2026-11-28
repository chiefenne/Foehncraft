# Install and update

Every release since 1.2.0 has ready-to-run programs for Linux, macOS and
Windows with `SHA256SUMS.txt` (README, "Download a release"); building from
source, as on this page, gives the newest `main`, MPI and the GPU backend.
For an install by an administrator under `/opt` that all users share, with
several versions side by side, see [Shared install on Linux](#shared-install-on-linux).
This page lists every tool to install on macOS,
Linux and Windows, with the commands, so that a fresh machine gets to a
working `foehn` in one pass. See also [Connect an assistant](agents.md) and
[Troubleshooting](troubleshooting.md).

Contents: [What you need](#what-you-need) ·
[macOS](#macos) · [Linux](#linux) · [Windows](#windows) ·
[Build](#build) · [Python, Gmsh and the example scripts](#python-gmsh-and-the-example-scripts) ·
[Optional: MPI and GPU](#optional-mpi-and-gpu) ·
[Shared install on Linux](#shared-install-on-linux) ·
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

These commands install Rust for your own user (rustup puts it into
`~/.cargo`). To install Foehncraft once for all users of a machine, as
root, see [Shared install on Linux](#shared-install-on-linux).

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
git checkout v1.2.0          # optional: the release instead of main
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

Then install what the scripts use (inside the virtual environment; a
`pip install` straight into the system Python stops on Debian 12 and newer,
Ubuntu 23.04 and newer and Fedora with `externally-managed-environment`):

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

## Shared install on Linux

For a machine (or a group of machines that mount the same folder) where an
administrator installs Foehncraft once, as root, and every user runs it:
the same idea as a commercial code under `/opt/<vendor>/<release>/`.
Several versions stay installed side by side, and each user picks one.

```
/opt/Foehncraft/
  1.2.0/bin/foehn            the program
  1.2.0/venv/                Python with gmsh, numpy, matplotlib, h5py
  1.2.0/src/                 the source tree of v1.2.0 (examples/, tools/), optional
  1.2.0/env.sh, env.csh      switch this version on (bash/zsh, csh/tcsh)
  1.3.0/...                  the next version, next to it
  current -> 1.2.0           the default version
  use.sh                     defines the command foehncraft-use (bash, zsh)
  modulefiles/foehncraft/    module files (Environment Modules, Lmod)
  rust/                      only for a build from source with --install-rust
```

Two things that go wrong when you follow the per-user steps as root:

- **rustup installs Rust for one user** (into `~/.cargo` of whoever runs
  it), so Rust installed as root is not found by other users. With the
  release programs no Rust is needed; for a build from source see
  [below](#build-from-source-into-opt).
- **`pip install` into the system Python stops with
  `error: externally-managed-environment`** on Debian 12 and 13, Ubuntu
  23.04 and newer, and Fedora (PEP 668). Use a virtual environment per
  version, or the distribution's packages, as below. Do not use
  `--break-system-packages`: it can break tools of the system itself.

### With the install script

[`tools/install-shared.sh`](../tools/install-shared.sh) does every step
below; `--dry-run` prints them without changing anything, `--help` lists
the options. As root:

```
curl -fsSLO https://raw.githubusercontent.com/chiefenne/Foehncraft/main/tools/install-shared.sh
bash install-shared.sh --dry-run                 # what it would do
bash install-shared.sh                           # the latest release
bash install-shared.sh --version 1.2.0 --with-source   # a given release, with examples/
```

What it does (latest release unless `--version`; prefix `/opt/Foehncraft`
unless `--prefix`):

1. downloads `foehncraft-<version>-linux-x86_64.tar.gz` and
   `SHA256SUMS.txt` from the release, checks the checksum, puts `foehn`
   into `<version>/bin/` (`--tarball FILE` installs a file you downloaded,
   for a machine without internet);
2. with `--with-source`, clones the tag into `<version>/src/` (the example
   scripts and cases are there, not in the release archive);
3. Python, `--python venv` (default): installs `python3-venv` and the
   libraries the Gmsh module loads (`libglu1-mesa libxcursor1 libxft2
   libxinerama1 libgomp1`, or the Fedora names) with apt or dnf, then
   makes `<version>/venv` and installs `gmsh numpy matplotlib h5py` there
   with pip; `--python system` installs the distribution packages instead
   (no pip at all); `--python none` skips Python;
4. writes `env.sh`, `env.csh`, the module file and `use.sh`, makes the new
   version the default (`current`; `--no-current` keeps the old default),
   makes everything readable for all users, and runs `foehn info`.

`--no-system-deps` leaves apt and dnf alone and only prints what is
missing. The release programs are built for x86-64 and need glibc 2.28 or
newer (RHEL, Rocky and Alma 8 and 9, Debian 11 to 13, Ubuntu 20.04 and
newer); v1.2.0 itself still needs glibc 2.39 (Debian 13, Ubuntu 24.04).
On older systems, or on ARM, build from source.

### Python for a shared install

**A virtual environment per version** (default, `--python venv`): the
system Python stays as it is, each Foehncraft version has the package
versions it was installed with, and users need no pip. By hand:

```
apt install python3-venv libglu1-mesa libxcursor1 libxft2 libxinerama1 libgomp1
python3 -m venv /opt/Foehncraft/1.2.0/venv
/opt/Foehncraft/1.2.0/venv/bin/pip install gmsh numpy matplotlib h5py
```

`env.sh` puts `venv/bin` on `PATH`, so `python3` is then that Python, with
the modules. Further packages that single scripts name (`scipy`, `pyvista`
...): `/opt/Foehncraft/1.2.0/venv/bin/pip install scipy`.

**The distribution's packages** (`--python system`), when pip is not
allowed on the machine; the modules then belong to the system `python3`:

```
# Debian and Ubuntu
apt install python3-numpy python3-matplotlib python3-h5py python3-gmsh
# Fedora
dnf install python3-numpy python3-matplotlib python3-h5py python3-gmsh
```

`apt show python3-gmsh` shows which Gmsh version the distribution has; it
is usually older than the pip one. A Miniconda or Anaconda Python of a user
does not matter for either way, as long as it is not first on `PATH` when
the scripts run (`conda deactivate`).

### Build from source into /opt

For an older Linux, ARM, MPI (`--features mpi`) or programs tuned for one
machine's CPU. The result has the same layout as above, so users switch
versions the same way. Complete for Debian 13 (Ubuntu the same), as root:

```
apt install build-essential cmake git curl python3 python3-venv \
            libglu1-mesa libxcursor1 libxft2 libxinerama1 libgomp1
# for --features mpi also: libopenmpi-dev openmpi-bin libclang-dev

curl -fsSLO https://raw.githubusercontent.com/chiefenne/Foehncraft/main/tools/install-shared.sh
bash install-shared.sh --version 1.2.0 --build --install-rust --native
```

- `--build` clones the tag into `<version>/src/` and runs
  `cargo install --locked --root /opt/Foehncraft/<version>` for `foehn`;
  with `--python venv` it also installs the
  [Python package](recipes/python.md#python-api) (`pip install ./python`,
  which builds the Rust module) into the venv. The build files go to a
  temporary folder that is removed at the end.
- `--install-rust` installs Rust system-wide into `/opt/Foehncraft/rust`
  (`RUSTUP_HOME` and `CARGO_HOME` there, readable for all), when no `cargo`
  is on `PATH`. A later `--build` uses it again; update it with
  `RUSTUP_HOME=/opt/Foehncraft/rust/rustup CARGO_HOME=/opt/Foehncraft/rust/cargo /opt/Foehncraft/rust/cargo/bin/rustup update stable`.
  Users never need it: the programs do not need Rust to run.
- `--native` optimises for the CPU of the build machine
  (`-C target-cpu=native`); the programs may then stop with "illegal
  instruction" on another CPU, so leave it out when the folder is shared by
  different machines.
- `--features mpi` builds the MPI version (`foehn run --ranks N`, see
  [MPI](#optional-mpi-and-gpu)).

Fedora, RHEL, Rocky and Alma: the same with
`dnf install gcc gcc-c++ make cmake git curl python3 mesa-libGLU libXcursor libXft libXinerama libgomp`
(and `openmpi-devel clang-devel` for MPI).

### Several versions side by side

Every version is a folder of its own; installing one never changes
another. Users switch in their shell:

```
source /opt/Foehncraft/use.sh       # once per shell, or at login, see below
foehncraft-use                      # list: 1.1.0, 1.2.0 (default) (active)
foehncraft-use 1.1.0                # this shell now runs 1.1.0
foehncraft-use default              # back to the default
```

With Environment Modules or Lmod:

```
module use /opt/Foehncraft/modulefiles
module avail foehncraft
module load foehncraft/1.1.0        # or: module load foehncraft (the default)
module switch foehncraft/1.2.0
```

csh and tcsh: `source /opt/Foehncraft/1.1.0/env.csh`. In a script or a
batch job, `source /opt/Foehncraft/1.1.0/env.sh` (or call
`/opt/Foehncraft/1.1.0/bin/foehn` directly) fixes the version for that job.

The administrator, as root:

```
ln -s /opt/Foehncraft/use.sh /etc/profile.d/foehncraft.sh   # default version on at every login
bash install-shared.sh --list                                # installed versions
bash install-shared.sh --activate 1.1.0                      # make 1.1.0 the default
```

A version is removed by deleting its folder and its module file
(`/opt/Foehncraft/1.1.0`, `/opt/Foehncraft/modulefiles/foehncraft/1.1.0`);
move the default to another version first if it pointed there.

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
| `error: externally-managed-environment` on `pip install` | pip into the system Python (PEP 668): use the virtual environment above, or for all users [a venv per version](#python-for-a-shared-install); never `--break-system-packages`. |
| `rustc` or `cargo` not found for a user after installing Rust as root | rustup installs for one user only: install as that user, or see [Build from source into /opt](#build-from-source-into-opt). |
| `version 'GLIBC_2.39' not found` when starting a release program | The Linux is older than that release's build: a newer release (glibc 2.28 from the next release on), or [build from source](#build-from-source-into-opt). |
