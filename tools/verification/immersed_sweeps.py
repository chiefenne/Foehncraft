"""The K/h sweeps of the penalised wall (L13 stage 1, D226; gallery 291).

Two sweeps in the mesh width h and the permeability K, with K = c h² (and
K = c h for the channel):

    channel   plane Poiseuille flow between two penalised walls that do not
              lie on cell faces (the channel of tests/immersed.rs, H = 0.74 m,
              16 to 256 cells across the box of height 1): the mass flow
              error and the position of the effective wall,
              (H − H_eff)/2 in cells with H_eff = (12 μ Q / (ρ G))^⅓
    sphere    Stokes flow through a simple cubic array of spheres
              (examples/cases/flow/immersed_stokes_sphere.toml, 16³ to 48³):
              the drag coefficient F/(6π μ a U) against the series of
              Hasimoto (1959) and Sangani & Acrivos (1982), 2.152

    python3 -m tools.verification.immersed_sweeps run FOEHN OUT.json   (about 1 h)
    python3 -m tools.verification.immersed_sweeps plot OUT.json PNG
"""
import json
import math
import os
import subprocess
import sys
import tempfile

H = 0.74
EXACT = 1.2 * H**3 / (12 * 0.1) * 0.1
C = 4 / 3 * math.pi * 0.2**3
X = C ** (1 / 3)
SERIES = 1 / (1 - 1.7601 * X + C - 1.5593 * C * C + 3.9799 * X**8 - 3.0734 * X**10)

CHANNEL = '''
case_format = "0.7"
[mesh.generate]
type = "hex"
dimension = 2
upper = [0.5, 1.0, 0.1]
cells = [4, {ny}, 1]
[fluid]
density = 1.0
viscosity = 0.1
pressure_values = "gauge"
[immersed]
permeability = {k:e}
[[immersed.body]]
name = "lower"
shape = "box"
lower = [-1.0, -1.0, -1.0]
upper = [2.0, 0.13, 1.0]
[[immersed.body]]
name = "upper"
shape = "box"
lower = [-1.0, 0.87, -1.0]
upper = [2.0, 2.0, 1.0]
[boundary.xmin]
type = "periodic"
partner = "xmax"
pressure_gradient = 1.2
[boundary.xmax]
type = "periodic"
[boundary.ymin]
type = "wall"
[boundary.ymax]
type = "wall"
[boundary.zmin]
type = "symmetry"
[boundary.zmax]
type = "symmetry"
[solver]
algorithm = "coupled"
tolerance = 1e-10
max_iterations = 3000
[output]
vtkhdf = false
wall_csv = false
'''


def solve(foehn, text, tmp):
    p = os.path.join(tmp, "case.toml")
    with open(p, "w") as f:
        f.write(text)
    out = os.path.join(tmp, "out")
    subprocess.run([foehn, "flow", p, "--out", out], capture_output=True, check=False)
    with open(os.path.join(out, "result.json")) as f:
        return json.load(f)["solves"][0]


def run(foehn, path):
    res = {"channel": [], "sphere": [], "series": SERIES}
    sphere = open(os.path.join(os.path.dirname(__file__), "../../examples/cases/flow/immersed_stokes_sphere.toml")).read()
    with tempfile.TemporaryDirectory() as tmp:
        for law, cs in (("h2", [1e-1, 3e-2, 1e-2, 1e-3, 1e-4]), ("h", [1e-3, 1e-4])):
            for c in cs:
                for ny in (16, 32, 64, 128, 256):
                    h = 1 / ny
                    k = c * h * h if law == "h2" else c * h
                    s = solve(foehn, CHANNEL.format(ny=ny, k=k), tmp)
                    m = s["periodic"][0]["mass_flow"]
                    heff = (12 * 0.1 * m / 0.1 / 1.2) ** (1 / 3)
                    res["channel"].append(dict(law=law, c=c, n=ny, converged=s["converged"],
                                               error=(m - EXACT) / EXACT, wall=(H - heff) / 2 / h))
                    print(res["channel"][-1], flush=True)
        for c in (1e-1, 3e-2, 1e-2, 1e-3):
            for n in (16, 24, 32, 48):
                h = 1 / n
                t = (sphere.replace("cells = [32, 32, 32]", f"cells = [{n}, {n}, {n}]")
                     .replace("permeability = 9.765625e-6", f"permeability = {c * h * h:e}")
                     .replace('derived = ["solid-fraction"]', "vtkhdf = false"))
                s = solve(foehn, t, tmp)
                u = s["periodic"][0]["mass_flow"]
                f = s["immersed"][0]["force"][0]
                kn = f / (6 * math.pi * 0.2 * u)
                res["sphere"].append(dict(c=c, n=n, converged=s["converged"], drag=kn,
                                          error=kn / SERIES - 1))
                print(res["sphere"][-1], flush=True)
    with open(path, "w") as f:
        json.dump(res, f, indent=1)


def plot(path, png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.ticker
    res = json.load(open(path))
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.6), dpi=150)
    colours = {1e-1: "#d9541e", 3e-2: "#e8a33d", 1e-2: "#1f5fa8", 1e-3: "#2a9d8f", 1e-4: "#7b4fa0"}
    for law, ls in (("h2", "-"), ("h", "--")):
        for c in sorted({r["c"] for r in res["channel"] if r["law"] == law}, reverse=True):
            rows = [r for r in res["channel"] if r["law"] == law and r["c"] == c and r["converged"]]
            if not rows:
                continue
            n = [r["n"] for r in rows]
            lab = f"K = {c:g} h²" if law == "h2" else f"K = {c:g} h (m)"
            axs[0].loglog(n, [abs(r["error"]) for r in rows], ls, marker="o", ms=3.5,
                          color=colours.get(c, "k"), label=lab)
            axs[1].semilogx(n, [r["wall"] for r in rows], ls, marker="o", ms=3.5,
                            color=colours.get(c, "k"), label=lab)
    n = [16, 256]
    axs[0].loglog(n, [0.2 * (16 / x) for x in n], ":", color="0.4", label="first order")
    axs[0].loglog(n, [0.2 * (16 / x) ** 2 for x in n], "-.", color="0.4", label="second order")
    axs[0].set_xlabel("cells across the box")
    axs[0].set_ylabel("|mass flow error|")
    axs[0].set_title("Penalised Poiseuille channel: mass flow", fontsize=10)
    axs[1].axhline(0, color="0.5", lw=0.8)
    axs[1].set_xlabel("cells across the box")
    axs[1].set_ylabel("wall position error (cells; > 0: wall inside the fluid)")
    axs[1].set_title("Penalised Poiseuille channel: effective wall", fontsize=10)
    for c in sorted({r["c"] for r in res["sphere"]}, reverse=True):
        rows = [r for r in res["sphere"] if r["c"] == c and r["converged"]]
        axs[2].semilogx([r["n"] for r in rows], [100 * r["error"] for r in rows], "-", marker="o",
                        ms=3.5, color=colours.get(c, "k"), label=f"K = {c:g} h²")
    axs[2].axhline(0, color="0.5", lw=0.8)
    axs[2].set_xlabel("cells per box edge (radius = 0.2 box)")
    axs[2].set_ylabel("drag coefficient error (%)")
    axs[2].set_title(f"Stokes sphere array, c = 0.0335: drag against {res['series']:.3f}", fontsize=10)
    for ax, ticks in zip(axs, ([16, 32, 64, 128, 256], [16, 32, 64, 128, 256], [16, 24, 32, 48])):
        ax.set_xticks(ticks, [str(t) for t in ticks])
        ax.minorticks_off() if ax is axs[2] else ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.grid(alpha=0.3, which="both")
        ax.legend(fontsize=7)
    fig.suptitle("Brinkman penalisation: sweeps in the mesh width h and the permeability K", fontsize=11)
    fig.tight_layout()
    fig.savefig(png, facecolor="white")
    print(f"wrote {png}")


if __name__ == "__main__":
    if len(sys.argv) == 4 and sys.argv[1] == "run":
        run(sys.argv[2], sys.argv[3])
    elif len(sys.argv) == 4 and sys.argv[1] == "plot":
        plot(sys.argv[2], sys.argv[3])
    else:
        sys.exit(__doc__)
