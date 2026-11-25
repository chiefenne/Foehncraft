"""Crops the README hero tiles from gallery figures (assets/readme/hero_*.jpg).

The tiles are cut from the gallery pictures, without axes and colour bars,
until 3D renders of these cases replace them one by one:

    python3 tools/plot/readme_hero.py

Run it again after a gallery figure below changes; the crop boxes are in
pixels of the figure as committed.
"""

from pathlib import Path

from PIL import Image

ROOT = Path(__file__).resolve().parents[2]
GALLERY = ROOT / "docs" / "figures" / "gallery"
OUT = ROOT / "assets" / "readme"

# (gallery figure, crop box left/top/right/bottom, output name, square pad)
TILES = [
    ("three_mixer_curved_rotor.png", (6, 53, 497, 544), "hero_curved_rotor.jpg", False),
    ("radiation_paths_3d.png", (112, 100, 815, 705), "hero_furnace_radiation.jpg", True),
    ("spray_two_way.png", (58, 34, 511, 487), "hero_spray_two_way.jpg", False),
    ("cht_heat_sink.png", (73, 40, 1218, 311), "hero_heat_sink.jpg", False),
]


def main() -> None:
    for name, box, out, square in TILES:
        tile = Image.open(GALLERY / name).convert("RGB").crop(box)
        if square:
            side = max(tile.size)
            canvas = Image.new("RGB", (side, side), "white")
            canvas.paste(tile, ((side - tile.width) // 2, (side - tile.height) // 2))
            tile = canvas
        tile.save(OUT / out, quality=90)
        print(f"{out}: {tile.width} x {tile.height}")


if __name__ == "__main__":
    main()
