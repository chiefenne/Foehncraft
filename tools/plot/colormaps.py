"""Colour maps of Foehncraft's pictures, the same as KoflCraft's and MurmelCraft's (D261).

Every field picture takes "Rainbow Desaturated"; Inferno is the second
map, for a second field next to the first (temperature next to the flow,
for example). A signed field whose range is centred on zero takes "Rainbow
Uniform" (Colin Ware): blue, green, a bright yellow at zero, orange, red,
magenta. The three projects share one definition: the control points
below, written out the same in MurmelCraft (``assets/colormaps.toml``) and
KoflCraft (``docs/colour-maps.md``, Rainbow Desaturated), and tested in each.

Source: the presets "Rainbow Desaturated" and "Rainbow Uniform" of
ParaView (``Remoting/Views/ColorMaps.json`` in github.com/Kitware/ParaView,
colour space RGB), Copyright (c) 2000 Kitware Inc., BSD 3-Clause licence
(``assets/licences/Rainbow-Desaturated-LICENSE.txt``). Every value of
Rainbow Desaturated is a whole number of 255ths, so its colours are exact;
Rainbow Uniform's are rounded to 8 bits, as in MurmelCraft.

    from tools.plot.colormaps import RAINBOW, UNIFORM  # registers both maps
    ax.pcolormesh(x, y, v, cmap=RAINBOW)
    ax.pcolormesh(x, y, w, cmap=UNIFORM, vmin=-m, vmax=m)

In ParaView (pvpython) the map is the built-in preset ``PARAVIEW_PRESET``,
or ``lut.RGBPoints = rgb_points(lo, hi)`` where the preset is missing.
Importing this module needs no matplotlib; the map is registered when
matplotlib is there.
"""

from __future__ import annotations

NAME = "rainbow_desaturated"
#: The name to pass as ``cmap``; ``RAINBOW + "_r"`` is the map reversed.
RAINBOW = NAME
#: The same map built into ParaView.
PARAVIEW_PRESET = "Rainbow Desaturated"
#: The second map (matplotlib's; ParaView calls it "Inferno (matplotlib)").
SECOND = "inferno"
#: The map of signed fields centred on zero, and its ParaView preset.
UNIFORM = "rainbow_uniform"
PARAVIEW_UNIFORM_PRESET = "Rainbow Uniform"

#: (position, colour) from the low end (0) to the high end (1), mixed in RGB.
POINTS = (
    (0.0, "#4747DB"),
    (0.143, "#00005C"),
    (0.285, "#00FFFF"),
    (0.429, "#008000"),
    (0.571, "#FFFF00"),
    (0.714, "#FF6100"),
    (0.857, "#6B0000"),
    (1.0, "#E04D4D"),
)

#: Rainbow Uniform: 43 colours at even steps from 0 to 1, mixed in RGB.
UNIFORM_POINTS = (
    "#0561FF",
    "#056CF7",
    "#0577F0",
    "#0582E8",
    "#058BDF",
    "#0595D5",
    "#059ECB",
    "#05A6BF",
    "#05AFB3",
    "#05B8A8",
    "#05C19A",
    "#05CA8C",
    "#05D47F",
    "#05DC6D",
    "#06E55C",
    "#04ED4B",
    "#46F327",
    "#7EF51C",
    "#A4F90C",
    "#C2FB09",
    "#E1FD06",
    "#FFFF03",
    "#FFF414",
    "#FFE826",
    "#FFDD37",
    "#FFD137",
    "#FFC437",
    "#FFB837",
    "#FFAC37",
    "#FF9F37",
    "#FF9337",
    "#FF8537",
    "#FF7737",
    "#FF6837",
    "#FE5536",
    "#FC4230",
    "#FD2636",
    "#F21E40",
    "#E6144A",
    "#DA0A54",
    "#CC0B5B",
    "#BD0C62",
    "#AE0D6A",
)


def rgb(colour: str) -> tuple[float, float, float]:
    """``#RRGGBB`` as three numbers from 0 to 1."""
    return tuple(int(colour[i : i + 2], 16) / 255 for i in (1, 3, 5))


def sample(t: float) -> tuple[float, float, float]:
    """The colour at ``t`` (0 low end, 1 high end, clamped), mixed in RGB."""
    t = min(max(t, 0.0), 1.0)
    for (t0, c0), (t1, c1) in zip(POINTS, POINTS[1:]):
        if t <= t1:
            f = (t - t0) / (t1 - t0)
            return tuple(a + f * (b - a) for a, b in zip(rgb(c0), rgb(c1)))
    return rgb(POINTS[-1][1])


def rgb_points(lo: float = 0.0, hi: float = 1.0) -> list[float]:
    """ParaView's ``RGBPoints`` (x, r, g, b, ...) of the map over [lo, hi]."""
    return [x for t, c in POINTS for x in (lo + t * (hi - lo), *rgb(c))]


def colormap(n: int = 1024):
    """The map as a matplotlib colour map of ``n`` entries."""
    from matplotlib.colors import LinearSegmentedColormap

    return LinearSegmentedColormap.from_list(NAME, [(t, rgb(c)) for t, c in POINTS], N=n)


def uniform_colormap(n: int = 1024):
    """Rainbow Uniform as a matplotlib colour map of ``n`` entries."""
    from matplotlib.colors import LinearSegmentedColormap

    return LinearSegmentedColormap.from_list(UNIFORM, [rgb(c) for c in UNIFORM_POINTS], N=n)


def register() -> str:
    """Registers both maps and their reverses with matplotlib (once); returns RAINBOW."""
    import matplotlib

    for name, make in ((NAME, colormap), (UNIFORM, uniform_colormap)):
        if name not in matplotlib.colormaps:
            cmap = make()
            matplotlib.colormaps.register(cmap)
            matplotlib.colormaps.register(cmap.reversed())
    return NAME


try:
    register()
except ImportError:  # pvpython without matplotlib: rgb_points still works
    pass
