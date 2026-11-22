"""Grid convergence of the turbulence models on the TMR flat plate and the
channel at Re_tau = 395, stage 6 part 6g, specification 9.3 and 9.4: the
observed order p, the extrapolated value and GCI_fine (Celik et al. 2008)
of the skin friction C_f(x = 0.97) and the drag C_D of the plate and of the
bulk velocity U_b+ of the channel, for every model on three meshes, and
the finest value against the references (White's correlation and the
Coles-Fernholz fit at the run's own momentum-thickness Reynolds number for
the plate, the DNS of Moser, Kim & Mansour for the channel). No number
computed by another solver is used (D201).

    python3 examples/studies/plate_channel_gci.py [OUT] [--models sst,zf,...]
        [--cases plate,channel] [--levels 0,1,2] [--jobs 4] [--report-only]

Each run is a separate `foehn flow` of the gallery case of the model on
one level, in OUT/<case>_<model>_<L> (default output/plate_channel_gci); a
finished run is not run again. The plate takes `verification.levels =
[L]` (the built-in mesher halves every interval in both directions L
times). The channel solution does not depend on x, so the channel is
refined in y only: its y nodes are those of level L (the base intervals
times 2^L, with the grading ratio that keeps the mesher's nested mapping,
`nested_ratio`), its two cells in x stay; on 768 cells across this is a
1 536-cell run instead of 12 288 with the same U_b+ (17.316 on both).

The channel runs converge to the scaled residual 10^-9
(`CHANNEL_TOLERANCE`) instead of the 10^-6 of the cases. On the fine
meshes the residuals pass 10^-6 while U_b+ is still moving: on 768 cells
across, SST stops 2.9 % low and EB 1.5 % high on 384 cells, which was the
"level-3 drop" of the first levels 1 to 3 run (2026-11-20). With 10^-9 the
stopped value is within 0.001 % of the converged one on every mesh. The
plate keeps the tolerance of its case. The models
that resolve the wall (SST low-re, k-ε-ζ-f compound, k-ε-ζ-f₀
compound-2023, EB compound) have a first-cell y+ about 1 and below. The
k-ε family runs with scalable wall functions, which take the wall cell at
y* = max(y*, 11.06): the refinement moves the first cell towards the wall
and the treatment keeps the result in the log layer, so the sweep checks
that this insensitivity holds (the GCI of a wall-function result, not of
the near-wall profile). A level whose run diverges or stops is reported.

Models (key, case files, treatment):

    sst   tmr_plate_sst / channel_sst              low-re
    keps  tmr_plate_keps / channel_keps            scalable
    rke   the same, realizable k-ε                 scalable
    zf    tmr_plate_zeta_f / channel_zeta_f        compound
    zf0   tmr_plate_zeta_f0 / channel_zeta_f0      compound-2023
    eb    tmr_plate_eb / channel_eb                compound

The finest k-ε plate (level 2 and above) and the k-ε-ζ-f plate and channel
of level 3 run with the relaxations 0.5 for momentum and turbulence: the
defaults diverge there (RUN014). Results are on the stage-6 verification
page: the plate on levels 0 to 3 (2026-11-19/20), the channel converged on
levels 0 to 3 (2026-11-20, PR #128, about 15 minutes on four cores):

    cargo build --release -p foehncraft-cli
    python3 examples/studies/plate_channel_gci.py output/gci_channel --cases channel --levels 0,1,2,3 --jobs 4
    python3 examples/studies/plate_channel_gci.py output/gci_channel --cases channel --levels 0,1,2 --report-only

Bands (specification 9.3, 9.4): plate GCI_fine of C_f ≤ 2 %, C_f within
±8 % of White and of Coles-Fernholz (the SST model itself sits about 6 %
below White, a finding, not a defect of the study); channel GCI_fine of
U_b+ ≤ 0.5 %, U_b+ within ±3 % of the DNS (17.55). Writes OUT/report.md
(also printed) and OUT/results.json. Needs the `foehn` executable
(FOEHN=target/release/foehn by default).
"""

from __future__ import annotations

import json
import math
import os
import subprocess
import sys
import tomllib
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CASES = ROOT / "examples/cases/flow"
FOEHN = os.environ.get("FOEHN", str(ROOT / "target/release/foehn"))

U_B = 17.55            # the bulk velocity of the channel cases (u_tau = 1 m/s)
UB_DNS = 17.55         # Moser, Kim & Mansour (1999), Re_tau 395
NU = 2.0e-7           # the plate cases: U = 1 m/s, nu = 2e-7 m²/s
RE_X = 0.970084 * 5e6  # Re_x of the probe (U L/nu = 5e6)
CF_WHITE = 0.455 / math.log(0.06 * RE_X) ** 2


def coles_fernholz(re_theta: float) -> float:
    """C_f = 2 [ln(Re_θ)/0.384 + 4.127]⁻², a fit to measured zero-pressure-
    gradient boundary layers (Nagib, Chauhan & Monkewitz 2007)."""
    return 2.0 / (math.log(re_theta) / 0.384 + 4.127) ** 2


def edge_profile(d: Path):
    """(u_e, θ) of the profile line `x097` of the run (U = 1 at the inlet):
    u_e the largest u of the line (the plate displaces the flow a little,
    +0.3 % at the edge of the layer), θ = ∫ (u/u_e)(1 − u/u_e) dy."""
    f = d / "line_x097.csv"
    if not f.exists():
        return None
    rows = [l.strip().split(",") for l in f.read_text().splitlines() if not l.startswith("#")]
    h = rows[0]
    y = [float(r[h.index("y")]) for r in rows[1:]]
    u = [float(r[h.index("u")]) for r in rows[1:]]
    ue = max(u)
    g = [(v / ue) * (1 - v / ue) for v in u]
    th = sum(0.5 * (g[i] + g[i + 1]) * (y[i + 1] - y[i]) for i in range(len(y) - 1))
    return ue, th



MODELS = {
    # key: (label, plate case, channel case, extra sets)
    "sst": ("k-ω SST, low-re", "tmr_plate_sst.toml", "channel_sst.toml", []),
    "keps": ("k-ε, scalable", "tmr_plate_keps.toml", "channel_keps.toml", []),
    "rke": ("realizable k-ε, scalable", "tmr_plate_keps.toml", "channel_keps.toml",
            ['turbulence.model="realizable-k-epsilon"']),
    "zf": ("k-ε-ζ-f, compound", "tmr_plate_zeta_f.toml", "channel_zeta_f.toml", []),
    "zf0": ("k-ε-ζ-f₀, compound-2023", "tmr_plate_zeta_f0.toml", "channel_zeta_f0.toml", []),
    "eb": ("elliptic-blending k-ε, compound", "tmr_plate_eb.toml", "channel_eb.toml", []),
}

# the first level that runs with the relaxations 0.5 (momentum and
# turbulence), per model and case
RELAXED_FROM = {("keps", "plate"): 2}

# iterations per level (fine meshes converge more slowly)
ITERATIONS = {0: 6000, 1: 12000, 2: 24000, 3: 40000}

# the plate runs of these models use the coupled algorithm (D178), within
# 3 000 iterations: ζ-f₀ converges on the 128 × 96 plate in 178 instead
# of 7 400 SIMPLE iterations, with C_f the same to 0.1 %
PLATE_COUPLED = {"sst", "zf", "zf0"}

# the channel: converged to this scaled residual with the coupled
# algorithm (D178), within these limits; the segregated SIMPLE of the case
# files needed up to 39 000 (ζ-f₀) and 91 000 (EB) iterations on 768
# cells across, the coupled algorithm 883 for ζ-f₀ (same U_b+ to five
# digits)
CHANNEL_TOLERANCE = 1e-9
CHANNEL_ITERATIONS = {0: 3000, 1: 3000, 2: 3000, 3: 3000}
# SST (1 536 cells) and EB (768 and 1 536) stop at this limit with a
# scaled residual near 10⁻⁷ (also with c_T 0.2), within 0.2 % of the U_b+
# of the converged SIMPLE runs of #128 (finding; the follow-up is the
# coupled algorithm, not longer runs)


def nested_ratio(ratio: float, n0: int, symmetric: bool, level: int) -> float:
    """The grading ratio of a mesh of n0·2^level intervals whose nodes are
    those of level `level` of the base mesh (n0 intervals, `ratio`): the
    mesher maps ξ with (e^{βξ} − 1)/(e^β − 1), β = ln(r) m0/(m0 − 1), m0
    the base intervals of one side (`Grading::nodes` in
    crates/foehncraft-mesh/src/generate.rs), so the same β on m0·2^level
    intervals is the ratio exp(β (m − 1)/m), m = m0·2^level."""
    m0 = n0 // 2 if symmetric else n0
    if level == 0 or abs(ratio - 1.0) < 1e-14 or m0 < 2:
        return ratio
    beta = math.log(ratio) * m0 / (m0 - 1)
    m = m0 << level
    return math.exp(beta * (m - 1) / m)


def channel_mesh(file: str, level: int) -> list[str]:
    """The sets of the channel mesh of `level`, refined in y only."""
    gen = tomllib.loads((CASES / file).read_text())["mesh"]["generate"]
    nx, ny, nz = gen["cells"]
    sets = [f"mesh.generate.cells=[{nx},{ny << level},{nz}]"]
    gy = gen.get("grading", {}).get("y")
    if gy is not None:
        r = nested_ratio(gy.get("ratio", 1.0), ny, gy.get("symmetric", False), level)
        sets.append(f"mesh.generate.grading.y.ratio={r!r}")
    return sets


def runs(models, cases, levels):
    """(directory name, case file, sets) of every run."""
    out = []
    for key in models:
        label, plate, channel, extra = MODELS[key]
        for case in cases:
            file = plate if case == "plate" else channel
            if file is None:
                continue
            for lv in levels:
                sets = ["verification.require_convergence=false", "output.vtkhdf=false"] + extra
                if case == "plate":
                    if key in PLATE_COUPLED:
                        sets += ["solver.max_iterations=3000", 'solver.algorithm="coupled"']
                    else:
                        sets += [f"solver.max_iterations={ITERATIONS[lv]}"]
                    sets += [f"verification.levels=[{lv}]"]
                else:
                    sets += [f"solver.max_iterations={CHANNEL_ITERATIONS[lv]}",
                             f"solver.tolerance={CHANNEL_TOLERANCE}",
                             'solver.algorithm="coupled"'] + channel_mesh(file, lv)
                if lv >= RELAXED_FROM.get((key, case), 99):
                    # the default relaxations (0.7, 0.8) diverge on the
                    # finest meshes (RUN014): k-ε in outer iteration 11 on
                    # the 128 × 96 plate (level 2), k-ε-ζ-f in iteration 51
                    # on the 512 × 384 plate and in iteration 1 182 on the
                    # level-3 channel
                    sets += ["solver.momentum.relaxation=0.5", "solver.turbulence.relaxation=0.5"]
                out.append((f"{case}_{key}_{lv}", file, sets))
    return out


def run_one(base: Path, item) -> int:
    name, file, sets = item
    d = base / name
    if (d / "flow.json").exists():
        return 0
    cmd = [FOEHN, "flow", str(CASES / file), "--out", str(d), "--threads", "1"]
    for s in sets:
        cmd += ["--set", s]
    d.mkdir(parents=True, exist_ok=True)
    with open(d.with_suffix(".log"), "w") as log:
        r = subprocess.run(cmd, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
    print(f"{name}: exit {r.returncode}", flush=True)
    return r.returncode


def values(d: Path, case: str) -> dict | None:
    if not (d / "flow.json").exists():
        return None
    s = json.loads((d / "flow.json").read_text())["solves"][-1]
    v = {"cells": s["cells"], "iterations": s["iterations"], "converged": s.get("converged")}
    if case == "plate":
        p = next(w for w in s["wall_probes"] if w["name"] == "x097")
        v["cf"] = p["skin_friction"]
        v["y_plus"] = p["y_plus"]
        e = edge_profile(d)
        v["re_theta"] = None if e is None else e[0] * e[1] / NU
        # C_f of the probe on the local edge velocity u_e, against the fit
        v["cf_cf"] = None if e is None else coles_fernholz(v["re_theta"])
        v["cf_local"] = None if e is None else v["cf"] / e[0] ** 2
        # C_D = F_x / (½ ρ U² · 2 m²) = F_x on the 1 m depth (ρ = U = 1)
        v["cd"] = next(b["force"][0] for b in s["boundaries"] if b["name"] == "plate")
    else:
        f = s["periodic"][0]["driving_force"]
        v["ub_plus"] = U_B / math.sqrt(f)
    return v


def celik(f1: float, f2: float, f3: float, r: float = 2.0) -> dict:
    """Observed order, extrapolated value and GCI_fine of three solutions,
    f1 the finest, constant refinement ratio r (Celik et al. 2008)."""
    e21, e32 = f2 - f1, f3 - f2
    if e21 == 0.0 or e32 == 0.0:
        return {"verdict": "unchanged", "p": None, "ext": f1, "gci": 0.0}
    q = e32 / e21
    if q < 0:
        # oscillatory: the band of the three solutions
        return {"verdict": "oscillatory", "p": None, "ext": None,
                "gci": 0.5 * (max(f1, f2, f3) - min(f1, f2, f3)) / abs(f1)}
    p = abs(math.log(abs(q))) / math.log(r)
    ext = (r ** p * f1 - f2) / (r ** p - 1.0)
    gci = 1.25 * abs(e21 / f1) / (r ** p - 1.0)
    return {"verdict": "monotonic" if q > 1.0 else "divergent", "p": p, "ext": ext, "gci": gci}


def pct(v, digits=2):
    return "–" if v is None else f"{100 * v:+.{digits}f} %".replace("-", "−")


def gci_txt(g):
    if g is None:
        return "–"
    return "< 0.01 %" if 100 * g < 0.01 else f"{100 * g:.2f} %"


def report(base: Path, models, cases, levels) -> str:
    res = {}
    lines = ["# Grid convergence of the turbulence models: flat plate and channel", "",
             f"Levels {', '.join(map(str, levels))}; GCI_fine after Celik et al. (2008), ratio 2.", ""]
    for case in cases:
        q = ("cf", "cd") if case == "plate" else ("ub_plus",)
        if case == "plate":
            lines += ["## TMR flat plate, Re_L = 5·10⁶", "",
                      f"References: C_f(0.97) White {CF_WHITE:.5e}; Coles-Fernholz at the Re_θ of the finest run of each model (Re_θ in the column).",
                      "", "| model | cells | y⁺ (finest) | C_f ·10³ | p | GCI_fine | C_f extrap. ·10³ | finest vs White | Re_θ (finest) | finest vs Coles-Fernholz | C_D ·10³ | p | GCI_fine | iterations, converged |",
                      "|---|---|---:|---|---:|---:|---:|---:|---:|---:|---|---:|---:|---|"]
        else:
            lines += ["", "## Channel, Re_τ = 395", "", f"Reference: U_b⁺ = {UB_DNS} (DNS).", "",
                      "| model | cells | U_b⁺ | p | GCI_fine | extrapolated | finest vs DNS | iterations, converged |",
                      "|---|---|---|---:|---:|---:|---:|---|"]
        for key in models:
            label, plate, channel, *_ = MODELS[key]
            if case == "channel" and channel is None:
                continue
            vs = [values(base / f"{case}_{key}_{lv}", case) for lv in levels]
            if any(v is None for v in vs):
                lines.append(f"| {label} | runs missing | | | | | | |")
                continue
            g = {k: celik(vs[-1][k], vs[-2][k], vs[-3][k]) if len(vs) >= 3 else None for k in q}
            res[f"{case}_{key}"] = {"values": vs, "gci": g}
            cells = " / ".join(str(v["cells"]) for v in vs)
            its = ", ".join(f"{v['iterations']}{'' if v['converged'] else '*'}" for v in vs)
            if case == "plate":
                cf = " / ".join(f"{1e3 * v['cf']:.4f}" for v in vs)
                cd = " / ".join(f"{1e3 * v['cd']:.4f}" for v in vs)
                gc, gd = g["cf"] or {}, g["cd"] or {}
                fin = vs[-1]["cf"]
                lines.append(
                    f"| {label} | {cells} | {vs[-1]['y_plus']:.2f} | {cf} | "
                    f"{'–' if gc.get('p') is None else format(gc['p'], '.2f')} ({gc.get('verdict', '')}) | {gci_txt(gc.get('gci'))} | "
                    f"{'–' if gc.get('ext') is None else format(1e3 * gc['ext'], '.4f')} | {pct(fin / CF_WHITE - 1)} | "
                    f"{'–' if vs[-1]['re_theta'] is None else format(vs[-1]['re_theta'], '.0f')} | "
                    f"{'–' if vs[-1]['cf_cf'] is None else pct(vs[-1]['cf_local'] / vs[-1]['cf_cf'] - 1)} | {cd} | "
                    f"{'–' if gd.get('p') is None else format(gd['p'], '.2f')} | {gci_txt(gd.get('gci'))} | {its} |")
            else:
                ub = " / ".join(f"{v['ub_plus']:.3f}" for v in vs)
                gu = g["ub_plus"] or {}
                lines.append(
                    f"| {label} | {cells} | {ub} | "
                    f"{'–' if gu.get('p') is None else format(gu['p'], '.2f')} ({gu.get('verdict', '')}) | {gci_txt(gu.get('gci'))} | "
                    f"{'–' if gu.get('ext') is None else format(gu['ext'], '.3f')} | {pct(vs[-1]['ub_plus'] / UB_DNS - 1)} | {its} |")
    lines += ["", "\\* the run stopped at its iteration limit before the tolerance (`converged` false)."]
    (base / "results.json").write_text(json.dumps(res, indent=1))
    return "\n".join(lines) + "\n"


def main(argv):
    args = [a for a in argv if not a.startswith("--")]
    opt = {}
    i = 0
    while i < len(argv):
        if argv[i].startswith("--"):
            if argv[i] == "--report-only":
                opt["report_only"] = True
            else:
                opt[argv[i][2:]] = argv[i + 1]
                args.remove(argv[i + 1]) if argv[i + 1] in args else None
                i += 1
        i += 1
    base = Path(args[0]) if args else ROOT / "output/plate_channel_gci"
    models = opt.get("models", ",".join(MODELS)).split(",")
    cases = opt.get("cases", "plate,channel").split(",")
    levels = [int(v) for v in opt.get("levels", "0,1,2").split(",")]
    base.mkdir(parents=True, exist_ok=True)
    if not opt.get("report_only"):
        todo = runs(models, cases, levels)
        # the finest meshes first, so that the long runs start early
        todo.sort(key=lambda t: -int(t[0].rsplit("_", 1)[1]))
        with ThreadPoolExecutor(int(opt.get("jobs", "4"))) as ex:
            codes = list(ex.map(lambda t: run_one(base, t), todo))
        if any(codes):
            print("some runs failed: see the .log files", flush=True)
    rep = report(base, models, cases, levels)
    (base / "report.md").write_text(rep)
    print(rep)


if __name__ == "__main__":
    main(sys.argv[1:])
