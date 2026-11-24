"""Showcase renderings of the shell and tube heat exchanger (application 203)
with ParaView: the hero view of the whole unit, a lengthwise mid cut with
the temperature of air, water and steel (conjugate heat transfer), the air
path around the baffles, and a close-up section at the middle baffle.

    /Applications/ParaView-6.2.0.app/Contents/bin/pvpython \
        examples/applications/shell_tube_hx/render_images.py \
        output/hx_counterflow/fields.vtkhdf docs/figures/applications [assets/readme]

Writes hx_render_hero.jpg, hx_render_cut.jpg, hx_render_airpath.jpg and
hx_render_baffle.jpg (1800 pixels wide); with a third argument also the
README teaser heat_exchanger_renderings.webp (the hero and the air path
cross-fading in a loop, made by `python3 render_images.py --teaser hero.jpg
airpath.jpg out.webp`: needs a Pillow with WebP). Off-screen with
pvpython (`--force-offscreen-rendering` on a headless box).

Look. Rasterised physically based shading (metallic and roughness) lit by
a procedural studio environment map, rendered at twice the size with a
transparent background and composited onto a white or dark studio
backdrop. OSPRay path tracing was tried first: in ParaView 6.2 on Apple
silicon its pathtracer loses the red channel (white comes out cyan), so it
is not used. Colour map "Rainbow Desaturated" everywhere; each range is the
2nd to 98th percentile of what the picture shows, so the variation inside
each region is visible; the water has its own temperature range.

Zones (`zone` cell data): 0 housing, 1 water (tube side), 2 air (shell
side), 3 tubesheets, 4 tubes, 5 baffles. Lengths in mm, the axis is y,
z is up.
"""
import os
import sys

import numpy as np
from PIL import Image, ImageFilter


def teaser(a, b, path, width=1200, hold=3.5, fade=0.7, steps=10):
    """README teaser: a and b (same size) alternate with a cross-fade, as an
    animated WebP that loops (full colour, unlike a GIF)."""
    h = round(a.size[1] * width / a.size[0])
    a, b = a.resize((width, h), Image.LANCZOS), b.resize((width, h), Image.LANCZOS)
    frames, durations = [], []
    for x, y in ((a, b), (b, a)):
        frames.append(x)
        durations.append(int(hold * 1000))
        for i in range(1, steps):
            frames.append(Image.blend(x, y, i / steps))
            durations.append(int(fade * 1000 / steps))
    frames[0].save(path, save_all=True, append_images=frames[1:], duration=durations,
                   loop=0, quality=86, method=6)
    print("wrote", path, "%dx%d, %d frames, %.0f kB" % (width, h, len(frames), os.path.getsize(path) / 1e3))


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "--teaser":
    # plain python3: the Pillow of pvpython has no WebP writer
    teaser(Image.open(sys.argv[2]).convert("RGB"), Image.open(sys.argv[3]).convert("RGB"), sys.argv[4])
    sys.exit(0)

from paraview import servermanager as sm
from paraview.simple import (
    Calculator, CellDatatoPointData, Clip, CreateTexture, CreateView, ExtractBlock,
    ExtractSurface, SurfaceNormals, GetColorTransferFunction, GetScalarBar,
    MergeBlocks, OpenDataFile, PointSource, Render, SaveScreenshot, Show, Slice,
    StreamTracerWithCustomSource, Threshold, Transform, Tube, Delete, Hide,
)
from vtkmodules.util.numpy_support import vtk_to_numpy

SRC = sys.argv[1]
OUT = sys.argv[2]
README = sys.argv[3] if len(sys.argv) > 3 else None
WIDTH = 1800
SS = 2  # supersampling factor
HERE = os.path.dirname(os.path.abspath(__file__))
SCRATCH = os.path.join(OUT, ".hx_render_tmp")
os.makedirs(SCRATCH, exist_ok=True)
CMAP = "Rainbow Desaturated"

# ---------------------------------------------------------------- data
reader = OpenDataFile(SRC)
internal = ExtractBlock(Input=reader, Selectors=["/Root/internal"])
grid = MergeBlocks(Input=internal)
mm = Transform(Input=grid)
mm.Transform.Scale = [1000.0, 1000.0, 1000.0]
calc = Calculator(Input=mm, AttributeType="Cell Data", ResultArrayName="T_C",
                  Function="temperature-273.15")
vel = Calculator(Input=calc, AttributeType="Cell Data", ResultArrayName="U",
                 Function="velocity_x*iHat+velocity_y*jHat+velocity_z*kHat")


def zones(lo, hi=None):
    """Cells of zones lo..hi as point data (each zone smoothed on its own)."""
    t = Threshold(Input=vel, Scalars=["CELLS", "zone"],
                  LowerThreshold=lo - 0.5, UpperThreshold=(hi if hi is not None else lo) + 0.5,
                  ThresholdMethod="Between")
    return CellDatatoPointData(Input=t, PassCellData=1)


housing, water, air = zones(0), zones(1), zones(2)
tubesheets, tubes, baffles = zones(3), zones(4), zones(5)
steel_bundle = zones(3, 5)


def fetch_array(src, name, assoc="POINTS"):
    d = sm.Fetch(src)
    if d.IsA("vtkMultiBlockDataSet") or d.IsA("vtkCompositeDataSet"):
        it = d.NewIterator(); it.InitTraversal(); vals = []
        while not it.IsDoneWithTraversal():
            b = it.GetCurrentDataObject()
            fd = b.GetPointData() if assoc == "POINTS" else b.GetCellData()
            a = fd.GetArray(name)
            if a is not None:
                vals.append(vtk_to_numpy(a))
            it.GoToNextItem()
        return np.concatenate(vals) if vals else np.zeros(0)
    fd = d.GetPointData() if assoc == "POINTS" else d.GetCellData()
    a = fd.GetArray(name)
    return vtk_to_numpy(a) if a is not None else np.zeros(0)


def prange(src, name, lo=2, hi=98, assoc="POINTS"):
    v = fetch_array(src, name, assoc)
    if v.ndim > 1:
        v = np.linalg.norm(v, axis=1)
    v = v[np.isfinite(v)]
    a, b = np.percentile(v, [lo, hi])
    return float(a), float(b)


# ---------------------------------------------------------------- look
def studio_env(path, bright=1.0):
    """Equirectangular studio environment for the image-based lighting.

    VTK takes the image with its zenith along +y and the azimuth u = 0 along
    -x, u = 1/4 along -z (checked with a mirror sphere); the studio is laid
    out in world coordinates (z up): a bright ceiling, a darker floor and
    soft boxes for the reflections in the steel."""
    W, H = 2048, 1024
    u = ((np.arange(W) + .5) / W)[None, :]
    th = ((np.arange(H) + .5) / H * np.pi)[:, None]
    a = 2 * np.pi * u
    d = np.stack(np.broadcast_arrays(-np.sin(th) * np.cos(a), np.cos(th), -np.sin(th) * np.sin(a)), -1)
    z = d[..., 2]
    img = (0.30 + 0.35 * np.clip(z, 0, 1) - 0.18 * np.clip(-z, 0, 1))[..., None] * np.ones(3)

    def box(direc, size, inten, col=(1, 1, 1)):
        n = np.array(direc, float)
        n /= np.linalg.norm(n)
        c = np.clip((d @ n), -1, 1)
        m = np.exp(-((np.arccos(c) / size) ** 6))
        img[:] += m[..., None] * inten * np.array(col)
    box((0.2, -0.3, 1.0), 0.55, 1.6)                 # ceiling panel
    box((1.0, -1.0, 0.7), 0.35, 2.2, (1, .97, .93))  # key, behind the camera
    box((-1.0, 0.6, 0.35), 0.25, 1.6, (.9, .95, 1))  # rim
    box((1.0, 1.0, 0.15), 0.15, 1.4)                 # strip low on the side
    box((-0.6, -1.0, 0.1), 0.20, 0.9)
    img = np.clip(img * bright, 0, 1) ** (1 / 2.2)
    Image.fromarray((img * 255).astype(np.uint8)).save(path)


ENV = os.path.join(SCRATCH, "studio_env.png")
studio_env(ENV)
env_tex = CreateTexture(ENV)


def view(size):
    v = CreateView("RenderView")
    v.ViewSize = [size[0] * SS, size[1] * SS]
    v.OrientationAxesVisibility = 0
    v.UseColorPaletteForBackground = 0
    v.BackgroundColorMode = "Skybox"
    v.BackgroundTexture = env_tex
    v.UseEnvironmentLighting = 1
    v.UseFXAA = 1
    v.KeyLightIntensity = 0.35
    v.UseAmbientOcclusion = 1
    v.AmbientSamples = 64
    return v


def steel(d, rough=0.32, aniso=0.0, colour=(0.78, 0.79, 0.81)):
    d.ColorArrayName = [None, ""]
    d.Interpolation = "PBR"
    d.Metallic = 1.0
    d.Roughness = rough
    d.Anisotropy = aniso
    d.DiffuseColor = list(colour)
    d.AmbientColor = list(colour)


def coloured(d, v, name, rng, title, rough=0.45, metallic=0.0, bar=None, lut_name=None):
    d.ColorArrayName = ["POINTS", name]
    lut = GetColorTransferFunction(lut_name or name)
    lut.ApplyPreset(CMAP, True)
    lut.RescaleTransferFunction(rng[0], rng[1])
    lut.AutomaticRescaleRangeMode = "Never"
    d.LookupTable = lut
    d.Interpolation = "PBR"
    d.Metallic = metallic
    d.Roughness = rough
    if bar is not None:
        sb = GetScalarBar(lut, v)
        sb.Title = title
        sb.ComponentTitle = ""
        sb.Visibility = 1
        sb.Orientation = "Horizontal"
        sb.WindowLocation = "Any Location"
        sb.Position = bar["pos"]
        sb.ScalarBarLength = bar.get("len", 0.22)
        sb.ScalarBarThickness = 14 * SS
        sb.TitleFontSize = 15 * SS
        sb.LabelFontSize = 13 * SS
        sb.TitleColor = bar["ink"]
        sb.LabelColor = bar["ink"]
        sb.TitleBold = 0
        ticks = nice_ticks(*rng)
        dec = max(0, -int(np.floor(np.log10(ticks[1] - ticks[0])))) if len(ticks) > 1 else 0
        sb.UseCustomLabels = 1
        sb.CustomLabels = ticks
        sb.LabelFormat = "{:.%df}" % dec
        sb.AddRangeLabels = 0
        sb.DrawTickMarks = 0
        sb.TitleFontFamily = "Arial"
        sb.LabelFontFamily = "Arial"
    return lut


def nice_ticks(a, b, n=5):
    """Round tick values inside [a, b]."""
    raw = (b - a) / n
    mag = 10 ** np.floor(np.log10(raw))
    step = min((m * mag for m in (1, 2, 2.5, 5, 10) if m * mag >= raw), default=10 * mag)
    t = np.arange(np.ceil(a / step) * step, b + 1e-9, step)
    return [float(x) for x in t]


def outline(v, cuts, shift, colour=(0.22, 0.24, 0.29), radius=0.22):
    """Thin dark tubes along the zone borders of a section, lifted a little
    towards the camera. All borders go into one actor: with one line or
    tube actor per zone the image-based lighting of the view fails and the
    whole picture comes out flat grey (ParaView 6.2, macOS)."""
    from paraview.simple import AppendGeometry, FeatureEdges, PassArrays
    edges = [PassArrays(Input=FeatureEdges(Input=c, BoundaryEdges=1, FeatureEdges=0,
                                           NonManifoldEdges=0, ManifoldEdges=0),
                        PointDataArrays=[], CellDataArrays=[], FieldDataArrays=[])
             for c in cuts]
    t = Transform(Input=AppendGeometry(Input=edges))
    t.Transform.Translate = list(shift)
    d = Show(Tube(Input=t, Radius=radius, NumberofSides=6), v)
    d.ColorArrayName = [None, ""]
    d.Interpolation = "PBR"
    d.Metallic = 0.0
    d.Roughness = 0.6
    d.AmbientColor = list(colour)
    d.DiffuseColor = list(colour)


def surface(src):
    return SurfaceNormals(Input=ExtractSurface(Input=src), FeatureAngle=40)


def backdrop(size, dark):
    W, H = size
    if not dark:
        y = np.linspace(0, 1, H)[:, None, None]
        top, bot = np.array([255, 255, 255]), np.array([238, 240, 244])
        a = top + (bot - top) * y ** 1.5
        return Image.fromarray(np.broadcast_to(a, (H, W, 3)).astype(np.uint8))
    yy, xx = np.mgrid[0:H, 0:W]
    r = np.hypot((xx - 0.5 * W) / (0.75 * W), (yy - 0.42 * H) / (0.75 * H))
    a = np.clip(1 - r, 0, 1) ** 1.4
    c0, c1 = np.array([10, 12, 17]), np.array([52, 58, 72])
    img = c0 + (c1 - c0) * a[..., None]
    return Image.fromarray(img.astype(np.uint8))


def save(v, name, size, dark, crop=True):
    """Colour from a render with the studio environment as background (the
    image-based lighting needs it), coverage from a second render with a
    transparent background; composited onto the backdrop."""
    res = [size[0] * SS, size[1] * SS]
    rgb_path = os.path.join(SCRATCH, name + "_rgb.png")
    tmp = os.path.join(SCRATCH, name + ".png")
    SaveScreenshot(rgb_path, v, ImageResolution=res)
    SaveScreenshot(tmp, v, ImageResolution=res, TransparentBackground=1)
    fg = Image.open(rgb_path).convert("RGB")
    fg.putalpha(Image.open(tmp).split()[3])
    bg = backdrop(fg.size, dark).convert("RGBA")
    # soft contact shadow under the object on the light backdrop
    if not dark:
        alpha = fg.split()[3]
        sh = Image.new("RGBA", fg.size, (40, 46, 60, 0))
        sh.putalpha(alpha.point(lambda a: int(a * 0.18)).filter(ImageFilter.GaussianBlur(18 * SS)))
        bg = Image.alpha_composite(bg, sh.transform(fg.size, Image.AFFINE, (1, 0, 0, 0, 1, -10 * SS)))
    img = Image.alpha_composite(bg, fg).convert("RGB")
    if crop:  # cut the empty margins (the colour bars count as content)
        x0, y0, x1, y1 = fg.split()[3].point(lambda a: 255 if a > 8 else 0).getbbox()
        pad = int(0.035 * img.size[0])
        img = img.crop((max(0, x0 - pad), max(0, y0 - pad), min(img.size[0], x1 + pad), min(img.size[1], y1 + pad)))
    w = WIDTH
    h = round(img.size[1] * w / img.size[0])
    img = img.resize((w, h), Image.LANCZOS)
    path = os.path.join(OUT, name + ".jpg")
    img.save(path, quality=92)
    print("wrote", path, "%dx%d" % (w, h))
    return img


# ---------------------------------------------------------------- streamlines
def air_streamlines(n_seeds, radius=13.0, x0=-85.0):
    """Streamlines from a disc of seeds just inside the air inlet nozzle."""
    seeds = PointSource(Center=[x0, 28.0, 23.0], Radius=radius, NumberOfPoints=n_seeds)
    st = StreamTracerWithCustomSource(Input=air, SeedSource=seeds)
    st.Vectors = ["POINTS", "U"]
    st.IntegrationDirection = "FORWARD"
    st.MaximumStreamlineLength = 3000.0
    st.IntegratorType = "Runge-Kutta 4-5"
    st.MaximumSteps = 20000
    return st


# ======================================================================
# 1. hero: the whole unit on a dark studio backdrop
size = (1600, 1000)
v = view(size)
# a quarter of the housing cut away (the quadrant x > 0, z > 0 facing the camera)
hs = Calculator(Input=surface(housing), ResultArrayName="q", Function="min(coordsX,coordsZ)")
hcut = Clip(Input=hs, ClipType="Scalar", Scalars=["POINTS", "q"], Value=0.0, Invert=1)
d = Show(hcut, v)
steel(d, rough=0.22, aniso=0.5)
bundle_s = surface(steel_bundle)
st_rng = prange(bundle_s, "T_C")
d = Show(bundle_s, v)
coloured(d, v, "T_C", st_rng, "Steel temperature (°C)", rough=0.38, metallic=0.25,
         lut_name="T_steel", bar=dict(pos=[0.05, 0.86], ink=[0.85, 0.87, 0.9]))
st = air_streamlines(400)
tb = Tube(Input=st, Radius=0.7, NumberofSides=10)
sp_rng = prange(st, "speed")
d = Show(tb, v)
coloured(d, v, "speed", sp_rng, "Air speed (m/s)", rough=0.5, lut_name="speed_air",
         bar=dict(pos=[0.33, 0.86], ink=[0.85, 0.87, 0.9]))
Render(v)
v.CameraPosition = [430, -260, 330]
v.CameraFocalPoint = [0, 150, 0]
v.CameraViewUp = [0, 0, 1]
v.CameraViewAngle = 26
hero_img = save(v, "hx_render_hero", size, dark=True, crop=False)
print("hero ranges: steel", st_rng, "speed", sp_rng)


# ======================================================================
# 2. lengthwise mid cut (plane x = 0): temperature of air, water and steel
size = (1600, 760)
v = view(size)
INK = [0.11, 0.12, 0.16]



def _cut(src, normal, origin):
    s = Slice(Input=src)
    s.SliceType = "Plane"
    s.SliceType.Origin = list(origin)
    s.SliceType.Normal = list(normal)
    s.Triangulatetheslice = 0
    return s


steel_all = zones(0, 0)
cut_air = _cut(air, (1, 0, 0), (0, 0, 0))
cut_water = _cut(water, (1, 0, 0), (0, 0, 0))
cut_steel = [_cut(z, (1, 0, 0), (0, 0, 0)) for z in (housing, tubesheets, tubes, baffles)]
a_rng = prange(cut_air, "T_C")
w_rng = prange(cut_water, "T_C")
s_vals = np.concatenate([fetch_array(c, "T_C") for c in cut_steel])
s_rng = tuple(float(x) for x in np.percentile(s_vals, [2, 98]))
d = Show(cut_air, v)
coloured(d, v, "T_C", a_rng, "Air temperature (°C)", rough=0.6, lut_name="T_air",
         bar=dict(pos=[0.06, 0.05], ink=INK))
d = Show(cut_water, v)
coloured(d, v, "T_C", w_rng, "Water temperature (°C)", rough=0.6, lut_name="T_water",
         bar=dict(pos=[0.39, 0.05], ink=INK))
for i, c in enumerate(cut_steel):
    d = Show(c, v)
    coloured(d, v, "T_C", s_rng, "Steel temperature (°C)", rough=0.5, lut_name="T_steel_cut",
             bar=dict(pos=[0.72, 0.05], ink=INK) if i == 0 else None)
outline(v, [cut_air, cut_water] + cut_steel, (0.3, 0, 0))
# the far half of the housing gives the cut its depth (the air nozzle at x = -100)
hb = Clip(Input=surface(housing), ClipType="Plane")
hb.ClipType.Origin = [0.5, 0, 0]
hb.ClipType.Normal = [1, 0, 0]
hb.Invert = 1
d = Show(hb, v)
steel(d, rough=0.25, aniso=0.4)
Render(v)
v.CameraPosition = [660, 150, 250]
v.CameraFocalPoint = [0, 150, 0]
v.CameraViewUp = [0, 0, 1]
v.CameraViewAngle = 27
save(v, "hx_render_cut", size, dark=False)
print("cut ranges: air", a_rng, "water", w_rng, "steel", s_rng)


# ======================================================================
# 3. the air path around the baffles: glass housing, steel bundle, dense streamlines
size = (1600, 1000)
v = view(size)
glass = surface(housing)
d = Show(glass, v)
d.ColorArrayName = [None, ""]
d.Interpolation = "PBR"
d.Metallic = 0.0
d.Roughness = 0.05
d.DiffuseColor = [0.85, 0.9, 1.0]
d.Opacity = 0.13
d = Show(surface(steel_bundle), v)
steel(d, rough=0.25, aniso=0.3)
st3 = air_streamlines(1400, radius=14.0, x0=-97.0)
tb3 = Tube(Input=st3, Radius=0.45, NumberofSides=8)
sp3 = prange(st3, "speed")
d = Show(tb3, v)
coloured(d, v, "speed", sp3, "Air speed (m/s)", rough=0.45, lut_name="speed_air3",
         bar=dict(pos=[0.05, 0.86], ink=[0.85, 0.87, 0.9]))
v.EnableRenderOnInteraction = 0
Render(v)
v.CameraPosition = [520, -170, 380]
v.CameraFocalPoint = [0, 150, 5]
v.CameraViewUp = [0, 0, 1]
v.CameraViewAngle = 25
airpath_img = save(v, "hx_render_airpath", size, dark=True, crop=False)
print("airpath range: speed", sp3)

# ======================================================================
# 4. close-up section through the middle baffle (y = 150 mm): temperature
#    of the baffle, tube walls and housing, water temperature in the tubes,
#    air speed in the baffle window
size = (1400, 1000)
v = view(size)
Y = 150.0
c_air = _cut(air, (0, 1, 0), (0, Y, 0))
c_water = _cut(water, (0, 1, 0), (0, Y, 0))
c_steel = [_cut(z, (0, 1, 0), (0, Y, 0)) for z in (housing, tubes, baffles)]
sp4 = prange(c_air, "speed")
w4 = prange(c_water, "T_C")
s4 = tuple(float(x) for x in np.percentile(np.concatenate([fetch_array(c, "T_C") for c in c_steel]), [2, 98]))
d = Show(c_air, v)
coloured(d, v, "speed", sp4, "Air speed (m/s)", rough=0.6, lut_name="speed_air4",
         bar=dict(pos=[0.10, 0.05], len=0.32, ink=INK))
d = Show(c_water, v)
d.ColorArrayName = [None, ""]
d.Interpolation = "PBR"
d.Roughness = 0.15
d.DiffuseColor = [0.55, 0.70, 0.86]  # water, plain: this picture is about the steel and the air
for i, c in enumerate(c_steel):
    d = Show(c, v)
    coloured(d, v, "T_C", s4, "Steel temperature (°C)", rough=0.5, lut_name="T_steel4",
             bar=dict(pos=[0.58, 0.05], len=0.32, ink=INK) if i == 0 else None)
outline(v, [c_air, c_water] + c_steel, (0, -0.3, 0))
# downstream of the cut: the housing and the bundle in steel for depth
for src in (housing, steel_bundle):
    cl = Clip(Input=surface(src), ClipType="Plane")
    cl.ClipType.Origin = [0, Y + 0.2, 0]
    cl.ClipType.Normal = [0, 1, 0]
    cl.Invert = 0
    d = Show(cl, v)
    steel(d, rough=0.25, aniso=0.3)
Render(v)
v.CameraPosition = [120, -95, 95]
v.CameraFocalPoint = [0, Y, 0]
v.CameraViewUp = [0, 0, 1]
v.CameraViewAngle = 27
save(v, "hx_render_baffle", size, dark=False)
print("baffle ranges: speed", sp4, "water", w4, "steel", s4)

if README:
    # the WebP writer is missing from the Pillow of pvpython: plain python3 makes the teaser
    import subprocess
    cmd = ["python3", os.path.abspath(__file__), "--teaser", os.path.join(OUT, "hx_render_hero.jpg"),
           os.path.join(OUT, "hx_render_airpath.jpg"), os.path.join(README, "heat_exchanger_renderings.webp")]
    if subprocess.run(cmd).returncode:
        print("teaser failed; run:", " ".join(cmd))
