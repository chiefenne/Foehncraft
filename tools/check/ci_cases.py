"""Runs one set of the case catalogue the way CI does (R5: each case checks
its own results and exits non-zero on failure).

    python -m tools.check.ci_cases SET [FOEHN]   run the set; FOEHN defaults to target/release/foehn
    python -m tools.check.ci_cases --check       every flow case file is in exactly one Linux set
    python -m tools.check.ci_cases --list SET    print the commands and case files of the set

Sets (times are the Linux CI runner, main run of 2026-11-15):
  linux-1  mesh, scalar, poisson, the 3D cavities           about 2.2 min
  linux-2  cylinder Re 20, cavity Re 1000, Poiseuille       about 2.7 min
  linux-3  Taylor-Green, transient cavity, the rest, and the
           chemistry reactors                               about 2.2 min
  long     cylinder Re 100 (long and short channels, statistics) and the
           cylinder Re 20 on three meshes                   main, nightly and full-ci only
  quick    macOS and Windows: mesh, scalar, poisson and the quick flows

Run from the repository root. Each case writes output/<case>/.
"""

from __future__ import annotations

import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
F = "examples/cases/flow/"

FLOW = {
    "linux-1": [
        "cavity_re100", "cavity_re100_columns", "cavity_re100_energy", "cavity_re100_log_steps", "cavity_re100_schedule", "cavity_re1000_poly3d", "couette_hex",
        "couette_mixed", "couette_poly", "couette_tet", "channel_slip", "poiseuille_table",
        "couette_derived", "couette_viscosity", "couette_transpiration", "poiseuille_water",
        "pipe_gmsh_raw", "channel_heated", "channel_sst", "channel_keps", "channel_zeta_f",
        "channel_zeta_f0", "cavity_re1000_coupled", "channel_heated_keps", "channel_heated_pr0025",
        "channel_eb", "channel_rotating_sst", "square_duct_qcr", "pipe_turbulent_sst",
        "tjunction_mixing", "diffuser_sst", "diffuser_channel_profile",
        "ubend_sst", "ubend_channel_profile",
    ],
    "linux-2": [
        "cylinder_re20", "cavity_re1000", "taylor_green_cfl10", "poiseuille_tet",
        "poiseuille_hex", "poiseuille_poly", "tmr_plate_sst", "tmr_plate_keps",
        "tmr_plate_zeta_f0", "diffuser_sst_sharp", "ubend_sst_sharp", "ranks_channel_poly",
        "ranks_heated_channel_sst", "ranks_restart_channel", "ranks_run_channel",
        "ranks_run_transient_sst", "ranks_run_mrf_couette",
        "ranks_run_sliding_couette",
        "ranks_run_imprint_couette",
        "immersed_cylinder_re20", "immersed_topology_bend", "immersed_mms",
        "immersed_cut_cylinder_re20",
        "ranks_water_if97_channel",
    ],
    "linux-3": [
        "species_binary_diffusion", "species_coaxial_jet", "flame_methane_air_jl4",
        "taylor_green_error", "adjoint_residual_cavity", "cylinder_potential", "cavity_re100_budget", "taylor_green_bdf2",
        "taylor_green_cn", "poiseuille_mixed", "cavity_re100_transient", "taylor_green_space",
        "cavity_re100_poly3d", "cavity_re100_tet3d", "pipe_wedge", "poiseuille_expression",
        "taylor_green_expression", "channel_ramp", "scalar_front_bdf2", "scalar_oblique_step_mlp",
        "poiseuille_periodic", "channel_heated_periodic", "pipe_heated_periodic",
        "pipe_rough_keps", "backstep_zeta_f0_coupled", "backstep_validation",
        "pipe_heated_keps", "plate_heated_keps", "jet_pipe_profile", "taylor_couette_mrf",
        "taylor_couette_sliding", "taylor_couette_imprint", "couette_two_rotors",
        "axial_stage_uniform", "pot_solid_body", "pot_solid_body_viscous",
        "naca0012_su2_re1000", "backstep_outlet_backflow", "taylor_couette_su2_zones",
        "cavity_heated", "conduction_mms", "couette_viscous_heating", "cht_composite_wall",
        "cht_heated_wall", "cht_turbulent_channel", "cht_contact_transient",
        "cht_heated_wall_transient", "cht_orthotropic_patch", "cht_circuit_board",
        "cht_heated_wall_accelerated", "cht_heated_wall_periodic_accelerated",
        "cht_pulsating_channel_cycle", "cht_counterflow_two_fluids", "cht_two_fluids_absolute",
        "cht_two_fluids_start",
        "cavity_sst_operator_check",
        "buoyancy_box_at_rest", "buoyancy_slot", "buoyancy_channel_reduced",
        "buoyancy_channel_static", "buoyancy_open_box", "buoyancy_enclosure_plume",
        "vd_heated_duct", "vd_closed_box_heating", "vd_closed_box_side_heated", "vd_mms_steady",
        "mixed_convection_channel",
        "spray_duct", "spray_datums", "spray_two_way", "spray_cooling", "spray_dispersion",
        "spray_cooling_steady", "spray_evaporation", "spray_heptane_evaporation",
        "porous_darcy_block", "porous_honeycomb_channel", "porous_schumann", "porous_storage_charge",
        "motion_free_stream", "motion_piston_adiabatic", "motion_channel_ranks",
        "steam_heated_duct", "steam_superheater_channel",
    ],
    "long": [
        "cylinder_re100", "cylinder_re100_short", "cylinder_re100_statistics", "cylinder_re100_convective", "cylinder_re20_gci",
        "tmr_plate_zeta_f", "backstep_zeta_f0", "tmr_plate_eb", "backstep_eb", "impinging_jet",
        "cht_heat_sink_turbulent", "cylinder_steady_diagnosis", "vd_cavity_large_dt",
        "cht_pulsating_channel_skin", "cht_heated_wall_frozen", "cht_pulsating_channel_warmup",
        "natural_convection_cavity_ra1e5", "natural_convection_plate",
        "flame_methane_air_gri30",
        "immersed_stokes_sphere", "spray_saturation_steady", "spray_breakup",
        "spray_chamber_ecn",
        "motion_piston_cooled_walls", "spray_g_g1",
    ],
    "quick": [
        "cavity_re100_budget", "couette_hex", "couette_mixed", "couette_poly", "couette_tet",
        "cavity_re100", "cavity_re100_poly3d", "poiseuille_hex", "poiseuille_poly",
        "taylor_green_space", "channel_slip", "poiseuille_table", "pipe_wedge", "couette_derived",
        "couette_viscosity", "poiseuille_water", "pipe_gmsh_raw", "scalar_front_bdf2",
        "poiseuille_periodic", "taylor_couette_mrf", "taylor_couette_sliding",
        "taylor_couette_imprint",
    ],
}
LINUX = ("linux-1", "linux-2", "linux-3", "long")

# Files named *_injectors.toml are included by a case (D246), not run.
# Cases whose expected exit code is not 0: the steady cylinder that sheds is
# stopped by its diagnosis on purpose (D198, exit code 3 = stopped).
EXPECTED_EXIT = {F + "cylinder_steady_diagnosis.toml": 3}
# Extra arguments of a case in CI: the ECN Spray G run (1.2 ms, about 30
# minutes on 4 threads) runs its first 0.1 ms there; the gallery has the
# whole run.
ARGS = {F + "spray_g_g1.toml": ["--set", "time.end_time=1e-4"]}


def _files(pattern: str) -> list[str]:
    """Case files below the repository root, sorted as `ls` sorts them."""
    return sorted(
        str(p.relative_to(ROOT)).replace("\\", "/")
        for p in ROOT.glob(pattern)
        if not p.name.endswith(("_gpu.toml", "_injectors.toml"))
    )


def _base() -> list[tuple[str, str]]:
    mesh = sorted(
        _files("examples/cases/mesh/*_box.toml") + _files("examples/cases/mesh/*_2d.toml")
    ) + ["examples/cases/mesh/poly_partition.toml", "examples/cases/mesh/naca0012_su2.toml"]
    return (
        [("mesh", f) for f in mesh]
        + [("scalar", f) for f in _files("examples/cases/scalar/*.toml")]
        + [("poisson", f) for f in _files("examples/cases/poisson/*.toml")]
    )


def cases(name: str) -> list[tuple[str, str]]:
    """The (command, case file) pairs of a set."""
    if name not in FLOW:
        raise SystemExit(f"unknown set '{name}' (sets: {', '.join(FLOW)})")
    flow = [("flow", F + c + ".toml") for c in FLOW[name]]
    if name == "linux-3":
        # The 0D reactors of stage SF2 (each checks its own `[reactor.check]`).
        flow += [("chemistry reactor", f) for f in _files("examples/cases/chemistry/*.toml")]
    return _base() + flow if name in ("linux-1", "quick") else flow


def check() -> int:
    """A new flow case that is in no set would silently never run: refuse it."""
    assigned = [f for s in LINUX for c, f in cases(s) if c == "flow"]
    dupes = sorted({f for f in assigned if assigned.count(f) > 1})
    missing = sorted(set(_files(F + "*.toml")) - set(assigned))
    absent = sorted(f for f in set(assigned) | {f for _, f in cases("quick")} if not (ROOT / f).is_file())
    if dupes:
        print("in more than one set: " + " ".join(dupes), file=sys.stderr)
    if missing:
        print("in no set (add to tools/check/ci_cases.py): " + " ".join(missing), file=sys.stderr)
    if absent:
        print("listed but missing: " + " ".join(absent), file=sys.stderr)
    if dupes or missing or absent:
        return 1
    print("every flow case is in exactly one Linux set")
    return 0


def run(name: str, foehn: str) -> int:
    status = 0
    for cmd, f in cases(name):
        print(f"::group::{cmd} {f}", flush=True)
        start = time.monotonic()
        code = subprocess.run([foehn, *cmd.split(), f, *ARGS.get(f, [])], cwd=ROOT).returncode
        if code != EXPECTED_EXIT.get(f, 0):
            print(f"::error::{cmd} {f} failed (exit code {code})", flush=True)
            status = 1
        print("::endgroup::")
        print(f"{cmd} {f}: {int(time.monotonic() - start)} s", flush=True)
    return status


def main(argv: list[str]) -> int:
    if not argv or argv[0] in ("-h", "--help"):
        print(__doc__, file=sys.stderr)
        return 0 if argv else 2
    if argv[0] == "--check":
        return check()
    if argv[0] == "--list":
        for cmd, f in cases(argv[1]):
            print(cmd, f, *ARGS.get(f, []))
        return 0
    foehn = argv[1] if len(argv) > 1 else "target/release/foehn"
    for name in (foehn, foehn + ".exe"):  # Windows: foehn.exe
        if (ROOT / name).is_file() or Path(name).is_file():
            foehn = str((ROOT / name) if (ROOT / name).is_file() else Path(name).resolve())
            break
    return run(argv[0], foehn)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
