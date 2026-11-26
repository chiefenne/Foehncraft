# Tools

The scripts around the solver, all in Python (D61). Run them from the
repository root as modules, `python3 -m tools.<group>.<name>`; each one
prints its usage with `--help` or in its docstring.

| Group | Scripts | What for |
|---|---|---|
| `tools/agent/` | `stdio_host` | a minimal MCP host that keeps `foehn mcp` open for scripts and assistants, with a transcript (L18, L19) |
| `tools/check/` | `ci_cases`, `compare_runs`, `backend`, `c_free`, `kernels`, `unsafe_modules`, `vtkhdf` | the case catalogue as CI runs it, comparison of runs (D31), the backend lint (D60), the C-free core (D87), the kernel rules (R10), `unsafe` only in the named modules (D202), VTKHDF files in VTK's reader |
| `tools/docs/` | `indexes`, `images` | the by-number tables of the gallery, recipe and application indexes; image sizes of the docs |
| `tools/bench/` | `build_hypre`, `cloud`, `g1_cloud`, `stage2_cloud`, `stage3_cloud`, `stage5_cloud`, `stage6_cloud`, `radiation_cloud`, `chemistry_cloud` | hypre from source (D20); the drivers for a cloud GPU (Colab, Kaggle, any Linux machine with an NVIDIA GPU), with their shared steps in `cloud` |
| `tools/mesh/` | `gmsh_test_meshes`, `su2_test_meshes`, `cgns_test_meshes`, `cgns_ngon_fixture`, `fluent_msh_fixtures`, `vtk_test_meshes`, `multizone_fixtures`, `taylor_couette_zones`, `render_mesh` | the test meshes and fixtures of the mesh readers; headless pictures of meshes and fields (L21) |
| `tools/plot/` | `*_chart`, `*_fields` and the other figure scripts, `backstep_paraview` | the charts and field images of the gallery and the verification pages; `backstep_paraview` runs in ParaView's Python (`pvpython tools/plot/backstep_paraview.py`) |
| `tools/release/` | `changelog`, `verification_manual`, `verification_status` | the changelog fragments into `CHANGELOG.md`; the generated verification manual (D68); the verification status page (D199) |
| `tools/verification/` | `vd_mms` | the manufactured solution of the variable density (stage 7c) |
| `tools/notebooks/` | `*_colab.ipynb` | Google Colab notebooks that run the cloud drivers on a free T4 |

`tools/requirements-ci.txt` lists the Python packages of the CI checks
(VTK, h5py, NumPy). The PTX kernels are built by
`python3 gpu/ptx-kernels/build.py`. The package `foehncraft` itself lives in
[`python/`](../python/README.md).
