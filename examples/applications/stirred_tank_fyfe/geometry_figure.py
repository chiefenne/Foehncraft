"""Vertical section of the 2 L tank of Fyfe et al. (the geometry of
`mesh.py`, from the constants of that script): vessel with the chamfered
bottom read from the Supporting Information text, liquid level, shaft, hub,
the blade zone of the pitched-blade turbine, the baffle plates and the
rotor zone of the frozen rotor. A side view, no top view (CONTRIBUTING.md,
images of four-blade rotors).

    python3 examples/applications/stirred_tank_fyfe/geometry_figure.py docs/figures/applications/fyfe_geometry.png
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon

sys.path.insert(0, str(Path(__file__).resolve().parent))
import importlib.util

spec = importlib.util.spec_from_file_location("fyfe_mesh", Path(__file__).with_name("mesh.py"))
src = Path(__file__).with_name("mesh.py").read_text(encoding="utf-8")
# read the constants only (mesh.py imports gmsh at the top)
ns = {}
for line in src.splitlines():
    if line[:2] in ("R ", "H ", "R_") or line.startswith(("CHORD", "T_BLADE", "Z_IMP", "N_BLADES", "Z_BAFFLE", "BAFFLE_R0", "Z0_ZONE")):
        pass
exec("import math\n" + "\n".join(l.split("#")[0] for l in src.splitlines()
                                  if l.startswith(("R ", "H ", "R_FLAT", "R_TIP", "CHORD", "T_BLADE", "Z_IMP", "R_HUB", "R_SHAFT", "R_ZONE", "R_COLUMN", "BAFFLE_R0", "Z_BAFFLE"))), ns)

R, H = ns["R"], ns["H"]
mm = 1000.0


def main():
    out = Path(sys.argv[1] if len(sys.argv) > 1 else "fyfe_geometry.png")
    fig, ax = plt.subplots(figsize=(5.2, 6.4), dpi=200)
    prof = [(-R, H), (-R, ns["R_FLAT"] * 0 + 0.015), (-ns["R_FLAT"], 0.0), (ns["R_FLAT"], 0.0), (R, 0.015), (R, H)]
    ax.add_patch(Polygon([(x * mm, z * mm) for x, z in prof], closed=False, fill=True, fc="#dbe9f6", ec="#1f4e79", lw=1.6))
    ax.plot([-R * mm, R * mm], [H * mm, H * mm], color="#1f4e79", lw=1.0, ls="--")
    ax.text(0.0, H * mm + 3, "liquid level (symmetry plane)", ha="center", fontsize=7, color="#1f4e79")
    # rotor zone (frozen rotor)
    ax.add_patch(Rectangle((-ns["R_ZONE"] * mm, 20), 2 * ns["R_ZONE"] * mm, 30, fill=False, ec="#999999", lw=0.8, ls=":"))
    ax.text(-ns["R_ZONE"] * mm, 51.5, "rotor zone", fontsize=6.5, color="#777777")
    # shaft and hub
    ax.add_patch(Rectangle((-ns["R_SHAFT"] * mm, 38), 2 * ns["R_SHAFT"] * mm, H * mm - 38, fc="#666666", ec="none"))
    ax.add_patch(Rectangle((-ns["R_HUB"] * mm, 26), 2 * ns["R_HUB"] * mm, 12, fc="#444444", ec="none"))
    # blades: projected height 15 mm, from the hub to the tip, 25 mm above the bottom
    for s in (-1, 1):
        x0, x1 = s * ns["R_HUB"] * mm, s * ns["R_TIP"] * mm
        ax.add_patch(Polygon([(x0, 25), (x1, 25), (x1, 40), (x0, 40)], closed=True, fc="#d9822b", ec="#8a4b10", lw=0.8))
    # baffle plates (plane through two of the four)
    for s in (-1, 1):
        a, b = sorted((s * ns["BAFFLE_R0"] * mm, s * (ns["BAFFLE_R0"] + 0.012) * mm))
        ax.add_patch(Rectangle((a, ns["Z_BAFFLE"] * mm), b - a, H * mm - ns["Z_BAFFLE"] * mm, fc="#7f9f7f", ec="#3b5b3b", lw=0.8))
    ax.annotate("chamfered bottom: cone from r = 7.5 mm,\n15 mm high at the wall (SI text, 15 x 20 mm\nthen 40 x 40 mm, read as one cone)", xy=(-40, 8.5),
                xytext=(-64, -12), fontsize=6.5, arrowprops=dict(arrowstyle="-", lw=0.6), ha="left", va="top")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(fc="#7f9f7f", ec="#3b5b3b", label="baffle plates (4, two in the plane)"),
                       Patch(fc="#d9822b", ec="#8a4b10", label="blades (4, 45°, D = 45 mm)"),
                       Patch(fc="#444444", label="hub and shaft")], loc="lower right", fontsize=6.5, bbox_to_anchor=(1.0, 0.0), frameon=False)
    ax.set_xlim(-R * mm - 6, R * mm + 6)
    ax.set_ylim(-40, H * mm + 12)
    ax.set_aspect("equal")
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("height above the bottom (mm)")
    ax.set_title("2 L tank of Fyfe et al.: vertical section through two baffle plates", fontsize=8)
    fig.tight_layout()
    fig.savefig(out)
    print("wrote", out)


if __name__ == "__main__":
    main()
