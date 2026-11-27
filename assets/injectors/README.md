# Standard injector tips

Two injector tips for the pictures of spray devices (D343, step 28 page
`docs/requirements/step28-injector-devices.md`). They are drawn by our own
script `tools/data/injector_standard.py` from published sizes and are used
by `geometry = "standard-gasoline"` and `geometry = "standard-diesel"` of a
`[[spray.device]]`. They are a picture, not a wall: the flow does not see
them (D344).

| File | Tip | Holes (default pattern) |
|---|---|---|
| `standard-gasoline.stl` | body 6 mm, domed tip of radius 4 mm | 8 holes of 0.165 mm, 0.17 mm long, drilled 37° from the axis, exits on a circle of radius 1 mm, hole 1 on +x |
| `standard-diesel.stl` | body 4 mm, sac dome of radius 1.2 mm with a cone tangent to it at 50° | 7 holes of 0.11 mm, 0.8 mm long, 77.5° from the axis (155° included spray angle), exits on a circle of radius 0.65 mm, hole 1 on +x |

## Frame and units

- Millimetres (the solver scales them by 1e-3 and by the device's `scale`).
- The tip at the origin, the injector axis along +z into the chamber, the
  azimuth zero along +x (hole 1). A device places the file by its frame:
  +z to `axis`, +x to the device's e₁ (from `reference`, turned by
  `clocking`), the origin to `mount`.
- Binary STL. Each hole is a bore (a tube of the hole diameter and length
  along the hole axis) and its mouth is cut out of the dome, so the holes
  show as dark openings. The surface is not closed at the hole mouths on
  purpose.

A case that gives its own `[[spray.device.hole]]` or `[[spray.device.pattern]]`
replaces the default holes; the drawn bores stay those of the file. A file
of your own: `geometry = "stl"`, `file = "my_tip.stl"` (relative to the file
that names it), `scale` (1e-3 for millimetres) and, when it is drawn in
another frame, `stl_origin`, `stl_axis`, `stl_reference` in the file's own
coordinates.

## Sources

- Gasoline: the hole count, diameter (0.165 mm), length (0.16 to 0.18 mm,
  0.17 mm used) and drill angle (37°) are the published Engine Combustion
  Network Spray G nozzle numbers as quoted in D309 (ECN, Spray G pages). The
  plume angle 33° is a case value of the injector, not of the tip. The body
  diameter, dome radius and exit circle radius are our choice (the exit
  radius is the one of `examples/cases/flow/spray_g_g1_injectors.toml`).
- Diesel: our choice, typical of a passenger-car common-rail nozzle (seven
  holes, 155° included spray angle); no published nozzle is copied. The
  single axial hole of ECN Spray A (0.09 mm, D310) is the pattern
  `count = 1, polar = 0` on the same body.
- The measured nozzle geometries of the ECN (X-ray tomography) are not
  stored here: their pages state no licence for redistribution
  (`docs/verification/DATA_SOURCES.md`). Download one yourself and give it
  as `geometry = "stl"` with its frame keys.

## Licence

Our own work, under the MIT licence of the repository. To redraw (and to
check that the files match the script):

    python3 tools/data/injector_standard.py
    python3 tools/data/injector_standard.py --check
