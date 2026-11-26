- **Install guide for macOS, Linux and Windows.** `docs/install.md` lists
  every tool with the install commands per system (Homebrew; apt and dnf;
  winget and Chocolatey): Rust, the C and C++ compiler, CMake (needed for
  the HDF5 and METIS builds), Git, Python, the Gmsh pip package for the
  example mesh scripts with the Linux libraries it needs (`libGLU` and X11),
  and Open MPI with libclang for `--features mpi`. The README's install
  section gives the short version per system. A 39 MB Gmsh wheel committed
  by accident in #196 is removed from the repository root.
