"""An animation of a field series in ParaView: one picture per
`fields_<step>.vtkhdf` of a transient run (`[output] fields_every`), seen
from +z on the mid-depth plane, with the walls drawn black, a fixed colour
range and the simulated time and rotor revolutions. Made for the stirred
vessel of gallery 98 to 100 (vorticity magnitude); any cell field of the
series works. Run it with ParaView's Python from the repository root:

    pvpython tools/plot/rotor_movie.py output/sliding_turbulent
    pvpython tools/plot/rotor_movie.py output/sliding_turbulent \\
        --field vorticity_magnitude --range 0 40 --size 1600 1600
    pvpython tools/plot/rotor_movie.py output/three_mixer \\
        --field vorticity --component z --range=-60:60 --size 1600 800

(on macOS: /Applications/ParaView-6.2.0.app/Contents/bin/pvpython). It
writes DIR/movie/frame_0000.png, ... and prints the ffmpeg command that
turns them into an MP4 (`brew install ffmpeg` on macOS). The picture fits
the domain to its shape (`--size 1600 800` for a wide vessel), and the label
gives the revolutions of every rotor of the case (A / B / C). A field name
that is not in the series stops the script with the list of the cell fields
that are.

Options:
    --field NAME     cell field (default vorticity_magnitude; e.g. speed, k,
                     pressure, relative_velocity_magnitude); a vector by its
                     base name (vorticity, velocity) with --component
    --component C    of a vector field: magnitude (default), x, y or z;
                     Foehncraft writes vectors as NAME_x, NAME_y, NAME_z and
                     NAME_magnitude, so `--field vorticity --component z`
                     is the cell array vorticity_z (the signed vorticity of a
                     2D case: use a range about zero, `--range=-60:60`)
    --range LO HI    fixed colour range, also as one argument `--range=LO:HI`;
                     pvpython turns an argument `-60` into `--60` before the
                     script sees it, so give a negative bound as
                     `--range=-60:60` (or `--range " -60" 60`) (default 0 to
                     the 99th percentile of
                     the last frame)
    --size W H       pixels (default 1600 1600)
    --out DIR        where the frames go (default DIR/movie)
    --every N        every N-th file only (default 1)
    --preset NAME    ParaView colour map by name (default "Inferno (matplotlib)";
                     e.g. "Cool to Warm", "Viridis (matplotlib)", "Turbo",
                     "Jet", "Rainbow Uniform", "Black-Body Radiation");
                     a name this ParaView does not have falls back to inferno
    --isocolor R G B colour of the isolines, 0 to 1
                     (default 1 1 1; on a light map use e.g. 0 0 0)
    --wallcolor R G B  colour of the wall lines (default 0 0 0)
    --background R G B  background (default 1 1 1)
    --colorbar on|off  the colour bar (default on)
    --cells          flat cell colours instead of the smoothed point values
                     (the default: cell data to point data, smooth contours)
    --isolines N     also N black contour lines of the field
    --slice          cut the mesh at mid-depth (3D cases); the default draws
                     the cells seen from +z, which also shows the polyhedral
                     cells of an imprinted sliding interface that a cut drops
    --fps N          frame rate of the printed ffmpeg command (default 25)
"""
import argparse
import glob
import math
import os
import re

from paraview import servermanager
from paraview.simple import (
    CellDatatoPointData, ColorBy, Contour, CreateView, ExtractBlock, FeatureEdges,
    GetColorTransferFunction, GetScalarBar, MergeBlocks, OpenDataFile, Render, SaveScreenshot, Show, Slice, Text,
)
from vtkmodules.util.numpy_support import vtk_to_numpy
import numpy as np


# inferno (matplotlib), 9 points
INFERNO = (
    (0.0000, 0.001, 0.000, 0.014), (0.1250, 0.120, 0.047, 0.282), (0.2500, 0.342, 0.062, 0.429),
    (0.3750, 0.569, 0.151, 0.404), (0.5000, 0.735, 0.216, 0.330), (0.6250, 0.867, 0.318, 0.227),
    (0.7500, 0.958, 0.469, 0.067), (0.8750, 0.987, 0.691, 0.090), (1.0000, 0.988, 0.998, 0.645),
)


def opt(obj, name, value):
    """Set a property that may not exist in every ParaView version."""
    try:
        setattr(obj, name, value)
    except Exception:  # noqa: BLE001  (also paraview.NotSupportedException)
        pass


def series(d, every):
    """The files of the series, by step, and their step numbers."""
    files = []
    for f in glob.glob(os.path.join(d, "fields_*.vtkhdf")):
        m = re.search(r"fields_(\d+)\.vtkhdf$", f)
        if m:
            files.append((int(m.group(1)), f))
    files.sort()
    return files[::every]


def case_times(d):
    """Time step, start time and the |Omega| of every rotor (in the order of
    the case file) from DIR/case.resolved.toml."""
    dt, t0, omega = None, 0.0, []
    try:
        with open(os.path.join(d, "case.resolved.toml")) as f:
            text = f.read()
    except OSError:
        return dt, t0, omega
    section = ""
    for line in text.splitlines():
        s = line.strip()
        if s.startswith("["):
            section = s.strip("[]")
            continue
        m = re.match(r"(\w+)\s*=\s*([-+0-9.eE]+)\s*$", s)
        if not m:
            continue
        key, val = m.group(1), float(m.group(2))
        if section == "time" and key == "time_step":
            dt = val
        elif section == "time" and key == "start_time":
            t0 = val
        elif section in ("sliding", "mrf") and key == "omega":
            omega.append(abs(val))
        elif section in ("sliding", "mrf") and key == "rpm":
            omega.append(abs(val) * 2.0 * math.pi / 60.0)
    return dt, t0, omega


def colour_range(args):
    """LO HI from `--range LO HI` or `--range=LO:HI` (also LO,HI); the
    values may carry spaces (`" -60"`, which pvpython leaves alone)."""
    parts = args[0].replace(",", ":").split(":") if len(args) == 1 else args
    try:
        lo, hi = (float(x) for x in parts)
    except ValueError:
        raise SystemExit(f"--range: two numbers, LO HI or LO:HI, not {' '.join(args)} "
                         "(pvpython turns -60 into --60: write --range=-60:60)")
    return lo, hi


def block(reader, name):
    """The blocks of the node `name` of a series file (`internal`: the
    cells, `boundary`: the wall surfaces), else None. The writer puts both
    in the data assembly; ParaView versions differ in whether a selector
    goes to the assembly or to the hierarchy of the blocks, so try both."""
    for assembly, selector in (("Assembly", f"//{name}"), ("Hierarchy", f"/Root/{name}"),
                               ("Hierarchy", f"//{name}")):
        part = ExtractBlock(Input=reader)
        try:
            part.Assembly = assembly
        except Exception:  # noqa: BLE001  (ParaView without the property)
            pass
        part.Selectors = [selector]
        part.UpdatePipeline()
        if part.GetDataInformation().GetNumberOfCells() > 0:
            return part
    return None


def fetch(src):
    """The cell data of a source, merged to one grid."""
    merged = MergeBlocks(Input=src)
    merged.UpdatePipeline()
    return servermanager.Fetch(merged).GetCellData()


def resolve(cells, field, component):
    """The cell array to draw and the component index of a multi-component
    array (None for a scalar): `field` itself, or of a vector written as
    separate scalars the array `field`_`component`."""
    names = [cells.GetArrayName(i) for i in range(cells.GetNumberOfArrays())]
    if field in names:
        arr = cells.GetArray(field)
        if arr.GetNumberOfComponents() == 1:
            return field, None
        if component == "magnitude":
            return field, -1
        return field, "xyz".index(component)
    if f"{field}_{component}" in names:
        return f"{field}_{component}", None
    raise SystemExit(f"no cell field `{field}` in the series; cell fields: {', '.join(names)}")


def values(cells, name, comp):
    """The values of the array `name` (its component `comp`, -1: magnitude)."""
    v = vtk_to_numpy(cells.GetArray(name))
    if comp is None:
        return v
    return np.linalg.norm(v, axis=1) if comp < 0 else v[:, comp]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dir")
    ap.add_argument("--field", default="vorticity_magnitude")
    ap.add_argument("--component", choices=["magnitude", "x", "y", "z"], default="magnitude")
    ap.add_argument("--range", nargs="+", default=None, metavar="LO HI")
    ap.add_argument("--size", type=int, nargs=2, default=[1600, 1600])
    ap.add_argument("--out", default=None)
    ap.add_argument("--every", type=int, default=1)
    ap.add_argument("--preset", default="Inferno (matplotlib)")
    ap.add_argument("--cells", action="store_true")
    ap.add_argument("--isolines", type=int, default=0)
    ap.add_argument("--isocolor", type=float, nargs=3, default=[1.0, 1.0, 1.0])
    ap.add_argument("--wallcolor", type=float, nargs=3, default=[0.0, 0.0, 0.0])
    ap.add_argument("--background", type=float, nargs=3, default=[1.0, 1.0, 1.0])
    ap.add_argument("--colorbar", choices=["on", "off"], default="on")
    ap.add_argument("--slice", action="store_true")
    ap.add_argument("--fps", type=int, default=25)
    a = ap.parse_args()
    files = series(a.dir, max(a.every, 1))
    if not files:
        raise SystemExit(f"no fields_<step>.vtkhdf in {a.dir} (set [output] fields_every)")
    out = a.out or os.path.join(a.dir, "movie")
    os.makedirs(out, exist_ok=True)
    dt, t0, omega = case_times(a.dir)

    reader = OpenDataFile(files[-1][1])
    interior = block(reader, "internal")
    if interior is None:
        # a reader without the assembly names: take everything
        interior = reader
        interior.UpdatePipeline()
    walls = block(reader, "boundary")
    info = interior.GetDataInformation()
    x0, x1, y0, y1, z0, z1 = info.GetBounds()
    zc = 0.5 * (z0 + z1)
    cells = fetch(interior)
    field, comp = resolve(cells, a.field, a.component)
    if a.range:
        lo, hi = colour_range(a.range)
    else:
        v = values(cells, field, comp)
        if v.min() < 0.0:   # a signed field (a component): symmetric about zero
            hi = float(np.percentile(np.abs(v), 99.0))
            lo = -hi
        else:
            lo, hi = 0.0, float(np.percentile(v, 99.0))
    print(f"{len(files)} frames, {field} from {lo:.4g} to {hi:.4g}")

    shown = interior
    if a.slice:
        shown = Slice(Input=interior)
        shown.SliceType = "Plane"
        shown.SliceType.Origin = [0.0, 0.0, zc]
        shown.SliceType.Normal = [0.0, 0.0, 1.0]
    if not a.cells:
        shown = CellDatatoPointData(Input=shown)
    if walls is not None:
        lines = Slice(Input=walls)
        lines.SliceType = "Plane"
        lines.SliceType.Origin = [0.0, 0.0, zc]
        lines.SliceType.Normal = [0.0, 0.0, 1.0]
    else:
        # no wall blocks found: the outline of the mid-depth cut of the
        # cells (vessel, blades and hubs; also an interface whose two sides
        # do not share their nodes)
        cut = Slice(Input=interior)
        cut.SliceType = "Plane"
        cut.SliceType.Origin = [0.0, 0.0, zc]
        cut.SliceType.Normal = [0.0, 0.0, 1.0]
        lines = FeatureEdges(Input=cut)
        opt(lines, "BoundaryEdges", 1)
        opt(lines, "FeatureEdges", 0)
        opt(lines, "NonManifoldEdges", 0)
        opt(lines, "ManifoldEdges", 0)

    v = CreateView("RenderView")
    v.ViewSize = list(a.size)
    v.OrientationAxesVisibility = 0
    v.Background = list(a.background)
    opt(v, "UseColorPaletteForBackground", 0)
    opt(v, "BackgroundColorMode", "Single Color")

    d = Show(shown, v)
    where = "CELLS" if a.cells else "POINTS"
    if comp is None:
        ColorBy(d, (where, field))
    else:
        ColorBy(d, (where, field, "Magnitude" if comp < 0 else "XYZ"[comp]))
    opt(d, "Ambient", 1.0)
    opt(d, "Diffuse", 0.0)
    lut = GetColorTransferFunction(field)
    # the preset names differ between versions (ParaView 6 has no
    # "Inferno (matplotlib)"): try the names, else the table of inferno
    applied = False
    try:
        from paraview.modules.vtkRemotingViews import vtkSMTransferFunctionPresets
    except Exception:  # noqa: BLE001
        from paraview.modules.vtkPVServerManagerRendering import vtkSMTransferFunctionPresets
    presets = vtkSMTransferFunctionPresets.GetInstance()
    for name in dict.fromkeys([a.preset, "Inferno (matplotlib)", "Inferno", "inferno"]):
        try:
            if presets.HasPreset(name) and lut.ApplyPreset(name, True):
                applied = True
                break
        except Exception:  # noqa: BLE001
            pass
    if not applied:
        lut.RGBPoints = [x for t, r, g, b in INFERNO for x in (lo + t * (hi - lo), r, g, b)]
        opt(lut, "ColorSpace", "RGB")
    opt(lut, "AutomaticRescaleRangeMode", "Never")
    lut.RescaleTransferFunction(lo, hi)
    bar = GetScalarBar(lut, v)
    for name, value in (("Title", field.replace("_magnitude", "").replace("_", " ")), ("ComponentTitle", ""),
                        ("TitleColor", [0.0, 0.0, 0.0]), ("LabelColor", [0.0, 0.0, 0.0]),
                        ("TitleFontSize", 14), ("LabelFontSize", 13), ("HorizontalTitle", 1),
                        ("WindowLocation", "Any Location"), ("WindowLocation", "AnyLocation"),
                        ("Orientation", "Vertical"),
                        ("Position", [0.86, 0.25]), ("ScalarBarLength", 0.5),
                        ("RangeLabelFormat", "%-#.3g"), ("AddRangeLabels", 1)):
        opt(bar, name, value)
    d.SetScalarBarVisibility(v, a.colorbar == "on")
    opt(bar, "Visibility", 1 if a.colorbar == "on" else 0)

    if a.isolines > 0 and not a.cells:
        mid = Slice(Input=shown)
        mid.SliceType = "Plane"
        mid.SliceType.Origin = [0.0, 0.0, zc]
        mid.SliceType.Normal = [0.0, 0.0, 1.0]
        iso = Contour(Input=mid, ContourBy=["POINTS", field],
                      Isosurfaces=[lo + (hi - lo) * (k + 1) / (a.isolines + 1)
                                   for k in range(a.isolines)])
        ci = Show(iso, v)
        ColorBy(ci, None)
        ci.AmbientColor = list(a.isocolor)
        ci.DiffuseColor = list(a.isocolor)
        opt(ci, "LineWidth", 1.5)
        opt(ci, "Position", [0.0, 0.0, 0.5 * (z1 - z0) + 1e-3 * (x1 - x0)])
    w = Show(lines, v)
    try:
        ColorBy(w, None)
    except Exception:  # noqa: BLE001  ("invalid association string 'NONE'")
        pass
    w.AmbientColor = list(a.wallcolor)
    w.DiffuseColor = list(a.wallcolor)
    opt(w, "LineWidth", 2.5)
    if not a.slice:
        # the cells are seen from +z: lift the lines of the mid-depth cut
        # above the top faces, else the blades are hidden under them
        opt(w, "Position", [0.0, 0.0, 0.5 * (z1 - z0) + 1e-3 * (x1 - x0)])

    label = Text(Text="")
    t = Show(label, v)
    opt(t, "WindowLocation", "Upper Left Corner")
    opt(t, "WindowLocation", "UpperLeftCorner")
    opt(t, "FontSize", 22)
    opt(t, "Color", [0.0, 0.0, 0.0])

    # ParaView resets the camera at the first render of a session: render,
    # then place the camera
    Render(v)
    xc, yc = 0.5 * (x0 + x1), 0.5 * (y0 + y1)
    v.CameraParallelProjection = 1
    v.CameraPosition = [xc, yc, zc + 10.0 * max(x1 - x0, y1 - y0)]
    v.CameraFocalPoint = [xc, yc, zc]
    v.CameraViewUp = [0.0, 1.0, 0.0]
    # the parallel scale is half the visible height: fit the domain to the
    # shape of the picture with a margin (and room for the colour bar)
    pw, ph = a.size
    room = 0.82 if a.colorbar == "on" else 0.96   # width share of the domain
    scale = 0.5 * max((y1 - y0) / 0.86, (x1 - x0) * ph / (room * pw))
    shift = 2.0 * scale * pw / ph * (0.5 - 0.5 * room)  # domain left of the bar
    v.CameraParallelScale = scale
    v.CameraFocalPoint = [xc + shift, yc, zc]
    v.CameraPosition = [xc + shift, yc, zc + 10.0 * max(x1 - x0, y1 - y0)]

    for i, (step, f) in enumerate(files):
        reader.FileName = f
        if dt is not None:
            time = t0 + step * dt
            s = f"t = {time:.3f} s"
            turns = [time * w / (2.0 * math.pi) for w in omega]
            if len(set(omega)) == 1:
                s += f",  {turns[0]:.2f} revolutions"
            elif turns:
                s += ",  revolutions " + " / ".join(f"{n:.2f}" for n in turns)
        else:
            s = f"step {step}"
        label.Text = s
        Render(v)
        SaveScreenshot(os.path.join(out, f"frame_{i:04d}.png"), v, ImageResolution=list(a.size))
        if i % 20 == 0 or i == len(files) - 1:
            print(f"  frame {i + 1} of {len(files)} (step {step})")
    print(f"wrote {len(files)} frames to {os.path.abspath(out)}")
    print("make the video with:")
    print(f"  ffmpeg -framerate {a.fps} -i {os.path.join(out, 'frame_%04d.png')} "
          f"-c:v libx264 -pix_fmt yuv420p -crf 18 {os.path.join(out, 'movie.mp4')}")


if __name__ == "__main__":
    main()
