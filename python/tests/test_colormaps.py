"""The colour map shared by Foehncraft, KoflCraft and MurmelCraft (D261):
ParaView's "Rainbow Desaturated", eight RGB control points written out the
same in the three repositories. A changed point or position breaks the
match between the three projects' pictures here."""

import sys

import pytest
from conftest import ROOT

sys.path.insert(0, str(ROOT))
from tools.plot import colormaps as cm  # noqa: E402

# The same table as MurmelCraft's assets/colormaps.toml and KoflCraft's
# docs/colour-maps.md.
POINTS = [
    (0.0, (0x47, 0x47, 0xDB)),
    (0.143, (0x00, 0x00, 0x5C)),
    (0.285, (0x00, 0xFF, 0xFF)),
    (0.429, (0x00, 0x80, 0x00)),
    (0.571, (0xFF, 0xFF, 0x00)),
    (0.714, (0xFF, 0x61, 0x00)),
    (0.857, (0x6B, 0x00, 0x00)),
    (1.0, (0xE0, 0x4D, 0x4D)),
]


def test_rainbow_desaturated_is_paraviews():
    assert cm.RAINBOW == "rainbow_desaturated"
    assert cm.PARAVIEW_PRESET == "Rainbow Desaturated"
    assert [t for t, _ in cm.POINTS] == [t for t, _ in POINTS]
    for (t, c), (_, want) in zip(cm.POINTS, POINTS):
        assert cm.rgb(c) == tuple(v / 255 for v in want), t
        assert cm.sample(t) == pytest.approx(tuple(v / 255 for v in want), abs=1e-12)
    # halfway between two points is their mean, mixed in RGB
    assert cm.sample(0.357) == pytest.approx((0.0, (255 + 128) / 510, 0.5), abs=1e-9)
    assert cm.sample(-1.0) == cm.sample(0.0)
    assert cm.sample(2.0) == cm.sample(1.0)


def test_paraview_points_span_the_range():
    p = cm.rgb_points(10.0, 20.0)
    assert len(p) == 32
    assert p[0] == 10.0 and p[-4] == 20.0
    assert p[4] == pytest.approx(11.43)
    assert p[1:4] == list(cm.rgb("#4747DB"))


def test_matplotlib_map_matches_the_points():
    matplotlib = pytest.importorskip("matplotlib")
    m = matplotlib.colormaps[cm.RAINBOW]
    for t, want in POINTS:
        got = [round(255 * v) for v in m(t)[:3]]
        # a table of 1024 entries: within two levels at the kinks
        assert all(abs(a - b) <= 2 for a, b in zip(got, want)), (t, got, want)
    assert matplotlib.colormaps[cm.RAINBOW + "_r"](0.0)[:3] == pytest.approx(cm.rgb("#E04D4D"))
