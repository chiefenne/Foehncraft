# Step 28: injector and nozzle devices (D340 to D347)

Maintainer request of 2026-11-27: "for injectors or other spray
introductions i want to have the possibility to visualize the injector or
nozzle device (general abstraction). injectors have an axis and location as
they are mounted, and their nozzle holes have orientations relative to
that. we need to finally have this represented somehow in our spray
configurations. dont change any current existing specification. i'd like to
provide an STL of a "standard injector" for gasoline and diesel and give the
user the possibility to specify their specific stl."

Companion of the [spray specification](step28-sprays.md) (D240 to D260,
D300 to D309, D320 to D331) and of the [Spray A page](step28-spray-a.md)
(D310 to D319). This page is a plan: nothing here is built yet, and coding
starts only on the maintainer's go. The four open points were answered by
the maintainer on decision cards on 2026-11-27 (section 7); the text follows
the answers.

## 1. What exists today (read on main, bebdb83)

- `[[injector]]` (D246, `foehncraft-particles/src/config.rs`,
  `struct Injector`) is **one hole**: an absolute `position` (m, where the
  parcels start) and `direction` (the plume axis), plus the nozzle keys of
  D306 and D307 (`hole_diameter`, `area_contraction`, `hole_length`,
  `discharge_coefficient`, `inlet_loss`, `rate`, `flow_share`) and the
  sampled injection of D322 (`samples`, `sector`). The struct has
  `deny_unknown_fields`, and so has the include file (`struct SprayFile`:
  `[liquid]` and `[[injector]]` only). A new key in either is a schema
  change, which is why everything below is a new optional table or a new
  optional key, never a changed meaning.
- A multi-hole injector is written as one `[[injector]]` per hole, each with
  the hole's absolute position and direction computed by hand. The Spray G
  file `examples/cases/flow/spray_g_g1_injectors.toml` (D309) holds hole 1
  of 8 only: exit 1 mm from the axis and 0.1 mm into the chamber, plume 33°
  from the injector axis, on the plume plane φ = −11.25°, moved 2 µm normal
  to the sector plane so that the parcels start inside the domain. The
  injector axis itself, the mounting point and the other seven holes exist
  only in the comment and in `tools/plot/spray_g_metrics.py`, which turns
  each parcel by k·45° for the projected liquid volume of the 8 plumes.
- Spray A (`spray_a_distortion_injectors.toml`) is a single axial hole on a
  quarter domain (`flow_share = 0.25`).
- STL is read in three places with the same reader
  (`foehncraft_mesh::format::stl::read`, ASCII and binary): `[[datum]]`
  surfaces (SP1b), `[[immersed.body]] shape = "stl"` (L13 stage 2, PR #322;
  cut cells) and the surface checks of `foehncraft_mesh::surface` (D225:
  welded corners, winding numbers, signed distance). Both STL keys take
  `scale` (default 1, 1e-3 for millimetres). None of them can place a file:
  there is no translation or rotation key.
- Output: parcels go to `spray/parcels_<step>.vtkhdf`; `series.json`
  (format `foehncraft.series` 1, D268) lists the files per step for the
  MurmelCraft time player; readers ignore keys they do not know (D80).
- Nothing in the output shows where the injector is: a viewer sees the
  parcels appear in empty space.

## 2. The idea in one paragraph

A **device** is the injector or nozzle body as mounted: a right-handed
frame (mounting point, axis, and a reference direction for the azimuth)
and, optionally, a surface (a standard STL or the user's own). Its **holes**
are given in that frame (angle from the axis, azimuth about it, exit radius
and setback), so a user describes a Spray G injector as "8 holes, 33° from
the axis, 45° apart" instead of 8 hand-computed vectors. An `[[injector]]`
can name a device hole; its position and direction are then computed from
the frame. Old files with absolute `position` and `direction` stay exactly
as they are, give the same results to the bit, and may still name a device
for the picture only.

## 3. Decisions (proposed)

| Id | Decision |
|---|---|
| D340 | **Device frame.** New optional table `[[device]]`, allowed in the case's `[spray]` section (`[[spray.device]]`) and in the injector include file (a new optional key of `SprayFile`, beside `[liquid]` and `[[injector]]`). Keys: `name`; `mount` (m, the tip point of the injector in the case coordinates, the origin of the frame); `axis` (the injector axis pointing from the body into the chamber, any length, normalised); `reference` (a direction not along `axis`; its part normal to the axis is the zero of the azimuth, default: the global axis least aligned with `axis`, the same rule as `sector` of D322). The frame is e₃ = axis, e₁ = reference made normal to the axis, e₂ = e₃ × e₁; an azimuth φ turns from e₁ towards e₂. Optional `clocking` (degrees, default 0) turns the whole device about its axis, so one file serves several mounting positions. Refused: zero or non-finite vectors, a reference along the axis, a name used twice. |
| D341 | **Holes in the device frame.** `[[device.hole]]` with `name`, `polar` (degrees, angle of the hole axis from the device axis), `azimuth` (degrees), `radius` (m, distance of the hole exit centre from the axis, default 0), `setback` (m, along the axis from `mount` into the chamber, default 0; negative is inside the body) and optional `exit_azimuth` (degrees) when the exit lies on another plane than the hole axis (default: `azimuth`). Hole exit x = mount + setback e₃ + radius (cos φₑ e₁ + sin φₑ e₂); hole direction n = cos θ e₃ + sin θ (cos φ e₁ + sin φ e₂). A ring of equal holes is `[[device.pattern]]` with `count`, `polar`, `radius`, `setback`, `first_azimuth` (default 0) and optional `names` (default `1`…`count`), which expands to `count` holes at first_azimuth + k·360°/count. Optional per-hole `diameter` (m) only for the picture and the consistency check (the flow keys stay on the injector). |
| D342 | **An injector names a hole** *(maintainer, card 1, 2026-11-27)*. `[[injector]]` gets two optional keys, `device` and `hole`. With them, `position` and `direction` may be left out and are computed from D341; with them **and** `position`/`direction` the absolute values win and the hole is used for the picture and a check (warning SPR0xx when the position is more than 5 % of the exit radius plus 10 µm from the hole exit, or the direction more than 0.5° from the hole axis; the Spray G file's 2 µm sector nudge passes). Optional `offset` (m, a vector in the case coordinates added to the computed exit) carries such nudges in the device form. `hole = "all"` expands one injector into one per hole of the device, each with the same keys and named `<injector>/<hole>`; `mass_flow` and `rate` are then per hole, as today. Without `device` nothing changes: `position` and `direction` stay required and mean what they mean today. The parcel file's `injector` index and every keyed random number (D121) are those of the expanded list, in the order of the file then the holes, so a case written either way gives the same bits when the list is the same. |
| D343 | **The surface of the device** *(frame: maintainer, card 4; source: card 3, section 6)*. `geometry` = `"none"` (default: the frame and the holes only, drawn as glyphs), `"standard-gasoline"`, `"standard-diesel"` or `"stl"` with `file` (relative to the file that names it, ASCII or binary) and `scale` (default 1, as the other STL keys). The STL is given in its own device frame: **tip at the origin, injector axis along +z into the chamber, azimuth zero along +x**, so the solver places it by the rotation [e₁ e₂ e₃] and the translation to `mount`. A file drawn in another frame gives `stl_origin`, `stl_axis` and `stl_reference` (in the file's own coordinates) instead of being re-exported. The standard files live in `assets/injectors/` with their frame, units and holes in a README; `"standard-*"` also gives default holes (`[[device.pattern]]`) which a case may replace. The surface is checked with D225 (open edges, flipped triangles: reported, not repaired) and must contain no hole exit inside the solid (warning). |
| D344 | **Visualisation only by default; the device as a wall later** *(maintainer, card 2, 2026-11-27)*. The device never changes the mesh or the flow unless `role = "wall"` is set. `role = "wall"` (later stage ID4) hands the placed surface to the cut-cell layer of L13 as an `[[immersed.body]] shape = "stl"` (which then needs no placement keys of its own) or to the mesh generator; it requires cells of a tenth of the tip size near the tip and is refused while the mesh cannot resolve the tip (checked with the cut-cell volume fractions). |
| D345 | **Output.** At the first output of a run, `spray/devices.vtkhdf` (and, without HDF5, `spray/devices.vtu`) with the placed surfaces (triangles, cell fields `device`, `part`), one line per hole from the exit along its axis (length `glyph_length`, default 5 hole diameters or 1 mm; fields `device`, `hole`, `injector`, `polar`, `azimuth`) and the cone of each injector's `half_angle` as an outline. The file is written once (the device does not move in this version) and listed in `series.json` under a new top-level key `"static": {"devices": "spray/devices.vtkhdf"}`; readers that do not know the key ignore it (D80), so the format stays version 1. `result.json` gets `spray.devices[]` with the frame and every computed hole exit and direction, so an agent reads the geometry back without parsing TOML. A sector or quarter model writes the whole device (all holes) and marks the holes whose injectors are in the domain (`in_domain`). |
| D346 | **Scriptable (D77 to D87).** Tool `spray.devices` of the registry (D196) and the MCP server (L18): given a case or an injector file, returns the frame, the expanded holes (exit, direction, injector names) and writes the placed surface and glyphs without running the flow (preview of a mounting before meshing). Same numbers as the run, same code path. `foehn check` lists the devices and the expanded injectors. Every new key is in `docs/keys.md` and the JSON schema; no interactive-only path. |
| D347 | **Not in this version.** Moving or needle-lifting devices (needle lift as a time table could later feed `rate` and the visualisation; mesh motion of MS0 would be needed for a wall); devices in the steady spray mode are allowed for the picture and the hole computation (no change to D260); nozzle-internal flow (the D307 `nozzle_k`/`nozzle_epsilon` stay the coupling point); radial-ring or swirl nozzles other than "polar and azimuth of a straight hole" (a hollow-cone swirl atomiser is one axial hole with a `half_angle` today and stays so). |

## 4. Examples (as they would read)

Spray G, all 8 holes from the device; the current hole-1 file stays as it
is and keeps giving the same bits:

```toml
# spray_g_device_injectors.toml (new file; spray_g_g1_injectors.toml untouched)
[[device]]
name = "spray-g"
geometry = "standard-gasoline"
mount = [0.0, 0.0, 0.0]
axis = [1.0, 0.0, 0.0]
reference = [0.0, 1.0, 0.0]

[[device.pattern]]
count = 8
polar = 33.0
radius = 1.0e-3
setback = 1.0e-4
first_azimuth = -11.25

[[injector]]
name = "g"
device = "spray-g"
hole = "1"                 # the sector model: one hole, half in the domain
offset = [0.0, 0.0, 2e-6]  # into the sector, as in the hand-written file
half_angle = 12.5
hole_diameter = 1.65e-4
flow_share = 0.5
# ... rate, area_contraction, hole_length, discharge_coefficient as today
```

The existing file with only a picture added, nothing else changed:

```toml
[[device]]                   # new lines
name = "spray-g"
geometry = "standard-gasoline"
mount = [0.0, 0.0, 0.0]
axis = [1.0, 0.0, 0.0]
reference = [0.0, 1.0, 0.0]

[[injector]]
name = "hole1"
device = "spray-g"           # new line: picture and check only
hole = "1"                   # new line
position = [1.000000e-04, 9.81175461e-04, -1.93128751e-04]   # unchanged
direction = [0.83867057, 0.53417395, -0.10625380]            # unchanged
```

A user's own diesel injector, 7 holes, mounted inclined in a bore:

```toml
[[device]]
name = "my-diesel"
geometry = "stl"
file = "my_injector_tip.stl"
scale = 1e-3                 # drawn in mm
mount = [0.0, 0.0, 0.0105]
axis = [0.0, 0.17365, -0.98481]   # 10° from -z
clocking = 12.0

[[device.pattern]]
count = 7
polar = 77.5                 # 155° included angle
radius = 0.6e-3
setback = 0.2e-3

[[injector]]
name = "d"
device = "my-diesel"
hole = "all"                 # seven injectors d/1 ... d/7
# nozzle keys per hole as today
```

## 5. Stages and pull requests

Sizes are agent working time for the code; CI (about 25 min per PR, CI
speed-up #288) comes on top. Each stage gets its gallery entry, case files
and recipe (D18), in the same PR; gallery and recipe numbers are the PR
numbers.

| Stage | Content | Size |
|---|---|---|
| ID0 | This plan (docs only). | done with this PR |
| ID1 | D340 to D342 and D346: `[[device]]`, `[[device.hole]]`, `[[device.pattern]]`, injector `device`/`hole`/`offset`, `hole = "all"`, the check warning, `foehn check` listing, tool `spray.devices` (JSON only), keys.md and schema. Tests: every existing spray case parses and runs to the same bits (the regression set of the spray tests unchanged); the device form of Spray G hole 1 gives the hand-written position within 2 µm (without `offset`) and to round-off with `offset`, the direction to 1e-8; a pattern of 8 holes reproduces the k·45° turns of `spray_g_metrics.py`; frames with a reference along the axis, zero vectors and unknown holes refused. Case: `spray_g_device` (all 8 holes on a full 360° coarse mesh, short run). Gallery: hole glyphs on the parcels. Recipe: "describe a multi-hole injector by its holes". | 30 to 45 min |
| ID2 | D343 (frame of the STL, placement, `stl_*` keys, D225 check) and D345 output (`spray/devices.vtkhdf`/`.vtu`, `series.json` `static`, `result.json` `spray.devices`), the tool writing the placed surface; ParaView script and `tools/plot` script that draw device + parcels (scripts list the available fields when one is missing). Gallery with 1-3 images: Spray G device with the 8 plumes, Spray A device with its plume; images looked at before merging (no four-fold or symbol-like patterns: the gasoline device shows 8 holes, so the picture is taken obliquely with the plumes coloured by speed, not end-on). Recipe: "show your own injector STL". | 45 to 60 min |
| ID2m | MurmelCraft: read `static.devices` from `series.json` and draw the device under the parcels in the time player (separate PR in MurmelCraft, after ID2 is merged). | 30 min |
| ID3 | The two standard STLs in `assets/injectors/` with README (frame, units, holes, source and licence) and `geometry = "standard-gasoline"` / `"standard-diesel"` with their default patterns, drawn by our script `tools/data/injector_standard.py` from published sizes (card 3, section 6). | 30 to 45 min |
| ID4 | D344 `role = "wall"`: the device surface as an immersed STL body (cut cells, L13 stage 2) or into the mesh generator; tip resolution check; a case with the plume next to the tip. Only with its own go. | not sized yet |

Order: ID1, ID2, ID3 (ID3 can come before ID2: the files do not depend on
the output), ID2m, ID4 last and separately decided.

## 6. Where the standard STLs come from *(maintainer, card 3, 2026-11-27)*

The maintainer chose: **we generate both standard files from published
sizes**, and the maintainer may replace them later with their own drawings
without any change of keys. The script `tools/data/injector_standard.py`
(our own geometry, parametric, in the repository with its output) draws:

- **standard gasoline**: a multi-hole direct-injection tip, a cylindrical
  body with a domed tip and a ring of straight holes. Default pattern from
  the published Spray G nozzle numbers quoted in D309: 8 holes of 0.165 mm,
  orifice length 0.17 mm, drill angle 37° (the plume angle 33° stays a case
  value of the injector), exit radius about 1 mm, with the body diameter and
  dome taken from published drawings of such tips (each number on the
  README with its citation, `not yet read` until checked).
- **standard diesel**: a common-rail tip with a sac and a ring of holes.
  Default pattern from published numbers of a multi-hole diesel nozzle (for
  example the ECN Spray B three-hole nozzle or a typical seven-hole
  passenger-car nozzle; the source is chosen and cited when ID3 is built);
  the single axial hole of Spray A (0.09 mm) is a pattern of `count = 1,
  polar = 0` on the same body.

Both are exported in the frame of D343 (tip at the origin, +z into the
chamber, azimuth zero along +x, metres) with a README in
`assets/injectors/` giving the frame, units, holes, every published number
with its citation and the MIT licence of the repository (the files are our
own work). A user's own STL replaces them with `geometry = "stl"`.

Not used as shipped files: the ECN nozzle geometries of Spray G and Spray A
(X-ray tomography by Argonne, published on the ECN pages). The pages state
no licence for redistribution, so by the rule of `docs/verification/DATA_SOURCES.md` (the repository
holds quoted numbers with their citation, not the data sets) they are not
stored in the repository; a user who downloads one gives it as
`geometry = "stl"` with its own frame keys, and a recipe shows how. The
published hole numbers (count, 33° plume angle, hole diameter) are quoted
with their citation as today (D309). Clean-room rule D201: the device is
geometry and placement, no numbers from another code are involved; the
standard tips are our own drawings from published sizes, not copies
of a commercial code's injector library.

## 7. Answers of the maintainer (cards, 2026-11-27)

1. **How the holes reach the injection** (D342): injectors may name a
   device hole and get their position and direction from it; the old
   absolute form stays as is. *Chosen as recommended.*
2. **Device in the flow** (D344): picture only now; the wall (cut cells)
   later, with its own go. *Chosen as recommended.*
3. **Who makes the standard STLs** (section 6): we generate both from
   published sizes. *The maintainer chose this instead of the recommended
   "the maintainer provides them, a stand-in until then".*
4. **Frame of an STL** (D343): one fixed convention plus the optional keys
   `stl_origin`, `stl_axis`, `stl_reference` for other CAD frames. *Chosen as
   recommended.*

Coding of ID1 to ID3 starts on the maintainer's go.

## 8. Checks (to be met by ID1 and ID2)

- VD.1 Existing spray cases unchanged: same parsed injectors, same parcel
  bits, same `result.json` apart from the new optional `spray.devices` (absent
  when no device is given).
- VD.2 Spray G hole 1 from the pattern: position within 2 µm of the
  hand-written file, to round-off with `offset`; direction within 1e-8.
- VD.3 Pattern turns: hole k equals hole 1 turned by k·360°/count about the
  axis (to round-off), the same turns as `spray_g_metrics.py`.
- VD.4 Placement: a placed STL's tip lies at `mount` and its +z on `axis`
  (to round-off); `stl_*` keys give the same placed triangles as a file
  re-exported in the convention.
- VD.5 `hole = "all"` gives the same parcels as the hand-written list of the
  same holes in the same order.
- VD.6 The tool `spray.devices` and the run write the same numbers.
