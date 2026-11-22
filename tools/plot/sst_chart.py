"""k-ω SST 2003 (stage 6 part 6b, D130, D134): the skin friction of the
flat plate along the plate against White's correlation and the
Coles-Fernholz fit to measured layers, the velocity profiles of the plate at x = 0.97 and of the
channel at Re_τ = 395 in wall units against Spalding's law, and the
first-cell y⁺ study of the automatic wall treatment. Prints the numbers and
writes the gallery chart.

    python3 -m tools.plot.sst_chart docs/figures/gallery/sst.svg PLATE CHANNEL YP2 YP5 YP10 YP30 YP100

PLATE is the output directory of tmr_plate_sst.toml on four meshes (the
finest is the reference of the third panel),
CHANNEL that of channel_sst.toml, YP* the plate on its second mesh with the
automatic treatment and the wall grading that gives a first-cell y⁺ of
about 2, 5, 10, 30 and 100 at x = 0.97:

    target/release/foehn flow examples/cases/flow/tmr_plate_sst.toml --out PLATE \\
        --set 'verification.levels=[0, 1, 2, 3]' --set 'solver.algorithm="coupled"'
    target/release/foehn flow examples/cases/flow/channel_sst.toml --out CHANNEL
    target/release/foehn flow examples/cases/flow/tmr_plate_sst.toml --out YP30 \\
        --set 'turbulence.wall_treatment="automatic"' --set 'mesh.generate.grading.y.ratio=361.1'

(ratios 7983.6, 2859.1, 1298.9, 361.1 and 83.3 for y⁺ ≈ 2, 5, 10, 30, 100).
"""
import csv
import math
import os
import sys

from tools.plot.property_chart import panel

X_PROBE = 0.970084                     # the probe of the plate cases
RE_L, NU_PLATE = 5e6, 2e-7
KAPPA, B = 0.41, 5.0
UB_DNS = 17.55                          # Moser, Kim & Mansour (1999), Re_τ = 395


def table(d, name):
    with open(os.path.join(d, name)) as f:
        rows = [r for r in f if not r.startswith("#")]
    return [{k: float(v) for k, v in r.items()} for r in csv.DictReader(rows)]


def white(re_x):
    return 0.455 / math.log(0.06 * re_x) ** 2


def coles_fernholz(re_theta):
    """C_f = 2 [ln(Re_θ)/0.384 + 4.127]⁻² (Nagib, Chauhan & Monkewitz 2007)."""
    return 2.0 / (math.log(re_theta) / 0.384 + 4.127) ** 2


def spalding_y(u):
    k = KAPPA * u
    return u + math.exp(-KAPPA * B) * (math.exp(k) - 1 - k - k * k / 2 - k ** 3 / 6)


def spalding_u(y):
    lo, hi = 0.0, 60.0
    for _ in range(100):
        m = 0.5 * (lo + hi)
        lo, hi = (m, hi) if spalding_y(m) < y else (lo, m)
    return 0.5 * (lo + hi)


def at(rows, x):
    r = min(rows, key=lambda r: abs(r["x"] - x))
    return 2 * r["tau_x"], r["y_plus"]


def main():
    if len(sys.argv) < 9:
        sys.exit(__doc__)
    dst, plate, chan = sys.argv[1:4]
    sweep = sys.argv[4:9]
    wall = sorted((r for r in table(plate, "wall_plate.csv") if r["x"] > 0), key=lambda r: r["x"])
    cf, _ = at(wall, X_PROBE)
    # Re_θ(x) from the momentum integral dθ/dx = C_f/2 (U = 1) along the plate
    theta, th_x, prev = 0.0, [], None
    for r in wall:
        c = 2 * r["tau_x"]
        if prev is not None:
            theta += 0.5 * (c + prev[1]) / 2 * (r["x"] - prev[0])
        th_x.append((r["x"], c, theta / NU_PLATE))
        prev = (r["x"], c)
    cf_fit = [(x, coles_fernholz(rt)) for x, _, rt in th_x if rt > 1000]
    print(f"plate: C_f({X_PROBE}) = {cf:.5e} on the finest mesh")
    for x in (0.2, 0.5, 0.97, 1.5, 1.9):
        c, _ = at(wall, x)
        print(f"  x = {x}: C_f = {c:.5e}, White {white(RE_L * x):.5e} ({100 * (c / white(RE_L * x) - 1):+.2f} %)")
    ut = math.sqrt(cf / 2)
    prof = [(r["wall_distance"] * ut / NU_PLATE, r["u"] / ut) for r in table(plate, "line_x097.csv")]
    prof = [(y, u) for y, u in prof if 0.1 <= y <= 3000]
    dev = max(abs(u / spalding_u(y) - 1) for y, u in prof if 1 <= y <= 300)
    print(f"plate x = 0.97: largest deviation of u⁺ from Spalding's law for 1 ≤ y⁺ ≤ 300: {100 * dev:.2f} %")
    cw = table(chan, "wall_ymin.csv")
    tau = sum(r["tau_x"] for r in cw) / len(cw)
    nu_c = 1 / 395
    utc = math.sqrt(tau)
    cprof = [(r["wall_distance"] * utc / nu_c, r["u"] / utc) for r in table(chan, "line_profile.csv")]
    cprof = [(y, u) for y, u in cprof if y >= 0.1]
    ub = 17.55 / utc
    print(f"channel: u_τ = {utc:.5f} m/s, Re_τ = {utc / nu_c:.1f}, U_b⁺ = {ub:.3f} (DNS {UB_DNS}, "
          f"{100 * (ub / UB_DNS - 1):+.2f} %), U_c⁺ = {max(u for _, u in cprof):.3f}")
    pts = []
    for d in sweep:
        c, yp = at(table(d, "wall_plate.csv"), X_PROBE)
        pts.append((yp, 100 * (c / cf - 1)))
        print(f"automatic, first-cell y⁺ {yp:7.2f}: C_f = {c:.5e} ({pts[-1][1]:+.2f} % against the finest mesh with low-re)")
    w, h = 1020, 380
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}" '
           'font-family="Helvetica, Arial, sans-serif" font-size="12">',
           f'<rect width="{w}" height="{h}" fill="white"/>']
    lg = math.log10
    panel(out, 70, 50, 250, 260, "flat plate, skin friction", "log₁₀ Re_x", "C_f × 10³", (5.7, 7.0), (2.0, 4.5),
          [("dash", [(lg(RE_L * x), 1e3 * white(RE_L * x)) for x in [0.1 * i for i in range(1, 21)]],
            "#888888", "White (2006)"),
           ("line", [(lg(RE_L * r["x"]), 1e3 * 2 * r["tau_x"]) for r in wall if r["x"] > 0.1], "#1f6fb4",
            "SST, 256 x 192 cells"),
           ("dash", [(lg(RE_L * x), 1e3 * c) for x, c in cf_fit if 0.2 <= x <= 2.0][::8], "#d9822b",
            "Coles-Fernholz (fit)")],
          ([6, 6.5, 7], [2, 2.5, 3, 3.5, 4, 4.5]), legend=110)
    sp = [(lg(y), spalding_u(y)) for y in [10 ** (i / 20) for i in range(-20, 71)]]
    panel(out, 410, 50, 250, 260, "velocity in wall units", "log₁₀ y⁺", "u⁺", (-1, 3.5), (0, 30),
          [("dash", sp, "#888888", "Spalding (κ 0.41, B 5.0)"),
           ("line", [(lg(y), u) for y, u in prof], "#1f6fb4", "plate, x = 0.97"),
           ("line", [(lg(y), u) for y, u in cprof], "#d9822b", "channel, Re_τ 395")],
          ([-1, 0, 1, 2, 3], [0, 5, 10, 15, 20, 25, 30]), legend=12)
    panel(out, 750, 50, 250, 260, "automatic treatment, C_f(0.97)", "log₁₀ first-cell y⁺",
          "C_f / C_f,finest − 1 (%)", (0, 2.2), (-4, 4),
          [("dash", [(0, 2), (2.2, 2)], "#888888", "±2 % band"),
           ("dash", [(0, -2), (2.2, -2)], "#888888", ""),
           ("line", [(lg(y), e) for y, e in pts], "#1f6fb4", "SST automatic"),
           ("dot", [(lg(y), e) for y, e in pts], "#1f6fb4", "")],
          ([0, 0.5, 1, 1.5, 2], [-4, -2, 0, 2, 4]), legend=12)
    out.append("</svg>")
    with open(dst, "w") as f:
        f.write("\n".join(out) + "\n")
    print(f"wrote {dst}")


if __name__ == "__main__":
    main()
