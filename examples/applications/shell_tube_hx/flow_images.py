"""Flow images of the heat exchanger (application 203): the temperature and
the speed on a section through both air nozzles and the axis, the
temperature on the plane x = 0 through the water nozzles and two tubes,
and the temperature of the steel (the housing cut open at x = 0 around the
whole tube bundle, and the bundle alone from the other side); with a third
argument also the README teaser.

    pip install h5py pyvista pillow
    xvfb-run -a python3 examples/applications/shell_tube_hx/flow_images.py \
        output/hx_coflow/fields.vtkhdf docs/figures/applications [assets/readme]

Off-screen; on a headless Linux box run it under xvfb-run. The script reads
the `VTKHDF/internal` group of `fields.vtkhdf` (cell data `temperature`,
`speed`, `zone`); zone 1 is the water (tube side), zone 2 the air (shell
side), the others are steel.
"""
import sys

import h5py
import numpy as np
import pyvista as pv
from PIL import Image
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[3]))
from tools.plot.colormaps import RAINBOW  # noqa: E402

pv.OFF_SCREEN = True
src, out = sys.argv[1], sys.argv[2]
readme = sys.argv[3] if len(sys.argv) > 3 else None

f = h5py.File(src, "r")["VTKHDF/internal"]
pts = f["Points"][()] * 1000.0  # mm
conn, off = f["Connectivity"][()], f["Offsets"][()]
g = pv.UnstructuredGrid(np.insert(conn, off[:-1], np.diff(off)), f["Types"][()], pts)
for k in ["temperature", "speed"]:
    g.cell_data[k] = f["CellData/" + k][()]
g.cell_data["zone"] = f["CellData/zone"][()].astype(int)
zone = g.cell_data["zone"]
fluid = (zone == 1) | (zone == 2)
g.cell_data["T_C"] = g.cell_data["temperature"] - 273.15
spd = g.cell_data["speed"].copy()
spd[~fluid] = np.nan
g.cell_data["speed_fluid"] = spd

BG = "#eef0f3"
INK = "#1b1f2a"
bar = dict(color=INK, title_font_size=30, label_font_size=26, n_labels=6, fmt="%.0f",
           width=0.5, position_x=0.25, position_y=0.04, height=0.06)


def plotter():
    p = pv.Plotter(off_screen=True, window_size=(2400, 1100))
    p.set_background(BG)
    p.enable_anti_aliasing("ssaa")
    return p


def crop(path, pad=30):
    """Cut the empty margins of a screenshot."""
    im = Image.open(path).convert("RGB")
    a = np.asarray(im).astype(int)
    bg = np.array([int(BG[i:i + 2], 16) for i in (1, 3, 5)])
    ys, xs = np.nonzero(np.abs(a - bg).sum(axis=2) > 12)
    box = (max(xs.min() - pad, 0), max(ys.min() - pad, 0),
           min(xs.max() + pad, im.width), min(ys.max() + pad, im.height))
    im.crop(box).save(path)


def cam(n, c, up, d=900.0):
    return [tuple(np.array(c) + d * np.asarray(n)), tuple(c), tuple(up)]


# section A: through both shell nozzles and the axis of the exchanger
nA = np.array([46.0, 0.0, 200.0])
nA /= np.linalg.norm(nA)
uA = np.cross(nA, [0, 1, 0])
uA /= np.linalg.norm(uA)
sA = g.slice(normal=nA, origin=(0, 0, 0))
centre = (0.0, 150.0, 0.0)

p = plotter()
p.add_mesh(sA, scalars="T_C", cmap=RAINBOW, clim=(15, 300),
           scalar_bar_args={**bar, "title": "temperature (°C)"})
p.camera_position = cam(nA, centre, uA)
p.camera.zoom(1.55)
p.add_text("air in, 300 °C, 15 m/s", position=(830, 960), font_size=14, color=INK)
p.add_text("air out, 1.5 bar", position=(1720, 170), font_size=14, color=INK)
p.add_text("water in the end chambers and tubes, 20 °C in", position=(60, 1030), font_size=14, color=INK)
p.screenshot(f"{out}/hx_flow_temperature.png")
crop(f"{out}/hx_flow_temperature.png")

p = plotter()
p.add_mesh(sA, scalars="speed_fluid", cmap=RAINBOW, clim=(0, 25), nan_color="#9aa1ad",
           scalar_bar_args={**bar, "title": "speed (m/s)"})
p.camera_position = cam(nA, centre, uA)
p.camera.zoom(1.55)
p.screenshot(f"{out}/hx_flow_speed.png")
crop(f"{out}/hx_flow_speed.png")

# section B: x = 0 through the water nozzles and two tubes
sB = g.slice(normal=(1, 0, 0), origin=(0, 0, 0))
p = plotter()
p.add_mesh(sB, scalars="T_C", cmap=RAINBOW, clim=(15, 300),
           scalar_bar_args={**bar, "title": "temperature (°C)"})
p.camera_position = cam([1.0, 0, 0], (0, 150, 25), (0, 0, 1))
p.camera.zoom(1.5)
p.screenshot(f"{out}/hx_flow_tube_plane.png")
crop(f"{out}/hx_flow_tube_plane.png")

# the steel: zones 0 (housing), 3 (tubesheets), 4 (tubes), 5 (baffles)
STEEL = (20, 30)
sbar = {**bar, "title": "steel temperature (°C)"}


def zones(grid, which):
    return grid.extract_cells(np.where(np.isin(grid.cell_data["zone"], which))[0])


bundle = zones(g, [3, 4, 5]).extract_surface(algorithm="dataset_surface")
house = zones(g, [0]).clip(normal=(1, 0, 0), origin=(0, 0, 0), invert=True)
p = pv.Plotter(off_screen=True, window_size=(2400, 1200))
p.set_background(BG)
p.enable_anti_aliasing("ssaa")
p.add_mesh(house.extract_surface(algorithm="dataset_surface"), scalars="T_C", cmap="inferno",
           clim=STEEL, scalar_bar_args=sbar)
p.add_mesh(bundle, scalars="T_C", cmap="inferno", clim=STEEL, show_scalar_bar=False, specular=0.2)
p.camera_position = [(470, -190, 300), (0, 150, 25), (0, 0, 1)]
p.camera.zoom(1.15)
p.screenshot(f"{out}/hx_steel_cutaway.png")
crop(f"{out}/hx_steel_cutaway.png")

p = pv.Plotter(off_screen=True, window_size=(2400, 1200))
p.set_background(BG)
p.enable_anti_aliasing("ssaa")
p.add_mesh(bundle, scalars="T_C", cmap="inferno", clim=STEEL, scalar_bar_args=sbar, specular=0.2)
p.camera_position = [(-470, -190, -250), (0, 150, 25), (0, 0, 1)]
p.camera.zoom(1.2)
p.screenshot(f"{out}/hx_steel_bundle.png")
crop(f"{out}/hx_steel_bundle.png")

if readme:
    # the README teaser: the temperature section without the labels, 3:1
    p = pv.Plotter(off_screen=True, window_size=(2400, 800))
    p.set_background(BG)
    p.enable_anti_aliasing("ssaa")
    p.add_mesh(sA, scalars="T_C", cmap=RAINBOW, clim=(15, 300), show_scalar_bar=False)
    p.camera_position = cam(nA, centre, uA)
    p.camera.zoom(1.6)
    tmp = f"{readme}/heat_exchanger_temperature.png"
    p.screenshot(tmp)
    crop(tmp, pad=40)
    im = Image.open(tmp).convert("RGB")
    # pad to 3:1 with the background, then scale to 1800 x 600
    w = max(im.width, 3 * im.height)
    h = w // 3
    canvas = Image.new("RGB", (w, h), tuple(int(BG[i:i + 2], 16) for i in (1, 3, 5)))
    canvas.paste(im, ((w - im.width) // 2, (h - im.height) // 2))
    canvas.resize((1800, 600), Image.LANCZOS).save(f"{readme}/heat_exchanger_temperature.jpg", quality=88)
    import os
    os.remove(tmp)
print("ok", sA.n_cells, sB.n_cells)
