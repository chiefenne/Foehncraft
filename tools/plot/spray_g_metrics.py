"""ECN Spray G metrics of a sector run (D309): the liquid and the vapour
penetration along the injector axis, measured as the ECN measures them.

The case computes one eighth of the 8-hole injector: a sector of 22.5° with
the plume on its side plane `zmin` (half a plume) and the plane half way to
the next plume as `zmax`, both mirror planes. The 8 plumes follow by the
symmetry of the injector: every parcel is mirrored about the plume plane and
the pair is turned by k·45°, k = 0..7.

Liquid (ECN, Pickett, ECN5.9): the projected liquid volume PLV(x, s) =
∫ LVF dy along the line of sight, in mm³ of liquid per mm² of image; the
liquid penetration is the largest distance x from the nozzle along the
injector axis where PLV exceeds a threshold, 0.2e-3 or 2e-3 mm³/mm². The
line of sight is perpendicular to the plane of plumes 1 and 5 (a side view
with two plumes in profile); the parcels' liquid volume is binned on pixels
of `PIXEL` (0.1 mm) and smoothed by a Gaussian of `SMOOTH` (0.1 mm), both
our choices: without smoothing the penetration moves by less than 0.1 mm
between 0.05 and 0.2 mm pixels, and a wider blur reaches past the last
parcels (0.6 mm farther with 0.25 mm at 0.1 ms).

Vapour: the largest distance along the injector axis where the fuel vapour
mass fraction exceeds `VAPOUR_THRESHOLD` (0.1 %, the ECN model definition of
the vapour penetration, our choice for Spray G as for Spray A), from the cells
of the sector, which hold the vapour of all plumes by symmetry.

    from tools.plot.spray_g_metrics import liquid_penetration, vapour_penetration
"""
import math

import numpy as np

SECTOR = 22.5          # degrees: one eighth of the injector, half a plume pitch
PLUMES = 8
PIXEL = 0.1            # mm
SMOOTH = 0.1           # mm, standard deviation of the Gaussian smoothing
VAPOUR_THRESHOLD = 1e-3


def full_injector(x):
    """The 8 plumes from the parcels of the sector (x, y, z in m, the plume on
    the plane φ = −SECTOR/2): 16 copies of every parcel, as an array of shape
    (16 n, 3), in a frame with plume 1 in the plane z = 0 (φ = 0)."""
    phi0 = -0.5 * math.radians(SECTOR)
    r = np.hypot(x[:, 1], x[:, 2])
    phi = np.arctan2(x[:, 2], x[:, 1]) - phi0      # 0 on the plume plane
    out = []
    for sign in (1.0, -1.0):                        # the mirror image
        for k in range(PLUMES):
            a = sign * phi + k * 2.0 * math.pi / PLUMES
            out.append(np.column_stack([x[:, 0], r * np.cos(a), r * np.sin(a)]))
    return np.concatenate(out)


def plv_image(x, d, n, xmax=0.08, half_width=0.06):
    """The projected liquid volume (mm³/mm²) of the full injector on pixels
    of PIXEL mm: axis 0 along the injector axis (0..xmax m), axis 1 across
    (−half_width..half_width m). The line of sight is z; plumes 1 and 5 lie
    in the image plane."""
    vol = n * math.pi / 6.0 * (d * 1e3) ** 3        # mm³ per parcel
    xs = full_injector(x) * 1e3                     # mm
    v = np.tile(vol, 2 * PLUMES)                    # the sector holds half a plume, its mirror the other half
    nx = int(round(xmax * 1e3 / PIXEL))
    ny = int(round(2 * half_width * 1e3 / PIXEL))
    img, _, _ = np.histogram2d(
        xs[:, 0], xs[:, 1], bins=[nx, ny],
        range=[[0.0, xmax * 1e3], [-half_width * 1e3, half_width * 1e3]], weights=v,
    )
    img /= PIXEL * PIXEL
    return smooth(img, SMOOTH / PIXEL)


def smooth(img, sigma):
    """A Gaussian blur of `sigma` pixels (separable, zero outside)."""
    if sigma <= 0.0:
        return img
    h = int(math.ceil(4.0 * sigma))
    k = np.exp(-0.5 * (np.arange(-h, h + 1) / sigma) ** 2)
    k /= k.sum()
    out = np.apply_along_axis(lambda r: np.convolve(r, k, mode="same"), 0, img)
    return np.apply_along_axis(lambda r: np.convolve(r, k, mode="same"), 1, out)


def liquid_penetration(x, d, n, threshold=0.2e-3):
    """The largest distance along the injector axis (mm) where the projected
    liquid volume exceeds `threshold` mm³/mm²."""
    img = plv_image(x, d, n)
    rows = np.nonzero((img > threshold).any(axis=1))[0]
    return 0.0 if len(rows) == 0 else (rows[-1] + 1) * PIXEL


def vapour_penetration(centres, y_fuel, threshold=VAPOUR_THRESHOLD):
    """The largest distance along the injector axis (mm) of a cell whose
    fuel vapour mass fraction exceeds `threshold`."""
    hit = y_fuel > threshold
    return 0.0 if not hit.any() else float(centres[hit, 0].max() * 1e3)
