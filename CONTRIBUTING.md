# Contributing

The standing rules (MUST and SHOULD) and the exact checks to run before
pushing are in [AGENTS.md](AGENTS.md); lessons from earlier work, including
experiments not to repeat, are in [docs/pitfalls.md](docs/pitfalls.md).
This page explains the documentation conventions in detail.

## Documentation for every feature (D18)

A feature is not done without three pieces of documentation, in the same
form as in PyAeroFlow:

1. **A gallery example**: one picture with the commands that reproduce it,
   in [docs/gallery.md](docs/gallery.md).
2. **Case files**: one or more commented case files in `examples/cases/`,
   listed in the [case catalogue](docs/cases.md) with the mesh they need
   and where they are discussed.
3. **A recipe**: a short command sequence for the task, naming the
   case-file keys it uses, in [docs/recipes.md](docs/recipes.md).

### Where the pages go

Each of the three is an index page with a one-line summary per topic page;
the entries themselves live on topic pages, each with its own table of
contents:

| Index | Topic pages | Split by |
|---|---|---|
| [docs/gallery.md](docs/gallery.md) | `docs/gallery/*.md` | feature (benchmarks, meshes, scalar and pressure, flow, boundary conditions, post-processing, run control) |
| [docs/recipes.md](docs/recipes.md) | `docs/recipes/*.md` | task area (benchmarks, development, meshes, scalar and pressure, GPU, flow, boundary conditions, post-processing, run control) |
| [docs/cases.md](docs/cases.md) | `docs/cases/*.md` | physics, following the folders of `examples/cases/` |

When you add an entry:

- put it on the topic page that fits; start a new topic page when a new
  physics or topic begins (energy, turbulence, reacting flows, sprays,
  multiphase, ...) and add a row for it to the index page;
- gallery examples and recipes take the number of the pull request that
  adds them (gallery 54, recipe 54; a second entry of the same pull request
  is 54b). Numbers are global and never reused, so "gallery 12" or
  "recipe 17" keeps its meaning in other pages, the changelog and commit
  messages; gaps are fine. Open the pull request first (a draft will do)
  to know the number. Entries up to 39 keep their consecutive numbers;
- give each entry an `<a id="short-name"></a>` anchor right above its
  heading `## 54. Title` and a row in its page's table of contents; the
  anchor names are unique per section and stay valid when an entry moves
  to another page;
- leave the "by number" table and the number column of the page table on
  `docs/gallery.md` and `docs/recipes.md` alone: the workflow
  `docs-index.yml` regenerates them after every merge to main
  (`python -m tools.docs.indexes`, `--check` to compare). A new topic page
  needs its row in the page table, with the numbers column left empty;
- link to an entry on its topic page (`gallery/flow.md#cavity`), not to the
  index; links to an anchor of the index page still work and land on the
  entry's row there;
- in the case catalogue, write the "Where" column as links
  (`[gallery 12](../gallery/flow.md#cavity)`), and add new keys to the
  page's key table; the complete key list `docs/keys.md` is generated from
  the code (`foehn schema`).

Keep a topic page below a few hundred lines; when one grows past that,
split it into two topic pages the same way.

### Images carry their size

Write an image on a gallery page or in the README as
`<img src="..." alt="..." width="W" height="H">`, or write `![alt](src)`
and run `python -m tools.docs.images`, which adds the size of the file
(scaled to at most 800 pixels wide, 200 in a table row). Without a size an
image takes no room until it has loaded, so a link to an entry far down a
page lands thousands of pixels too high while the images above it arrive.
CI runs `python -m tools.docs.images --check`.

### Colour maps

Field pictures use the colour maps shared with KoflCraft and MurmelCraft
(D261): **Rainbow Desaturated** for every field, **Inferno** as the second
map (a second field beside the first, such as the temperature next to the
flow), and **Rainbow Uniform** (Colin Ware) for signed fields whose range
is centred on zero (zero is its bright yellow middle). In matplotlib
scripts import the maps from `tools/plot/colormaps.py`
(`from tools.plot.colormaps import RAINBOW, UNIFORM`, then `cmap=RAINBOW`
or `cmap=UNIFORM`); in ParaView use the presets `"Rainbow Desaturated"`
and `"Rainbow Uniform"`, or `colormaps.rgb_points(lo, hi)` where the first
is missing. Line colours
(streamlines over a field, series in a chart) may use other maps.

### Images must not resemble symbols

Never publish an image (gallery, README, docs, posts) whose pattern could
be misread as a symbol. Flow fields of rotating machinery take the
symmetry of the rotor, and a four-blade rotor in particular draws a
four-armed pinwheel in its vorticity and speed that can be mistaken for a
hateful symbol. So:

- show rotors with **three or five blades** in published images (the
  gallery vessels: `examples/meshes/mixer.py --blades 5 --baffles 3`,
  `examples/meshes/twin_mixer.py --right-blades 5`); four-blade meshes,
  cases and tests may stay, and their numbers may appear as text and in
  line plots, but not their flow fields;
- avoid other four-fold, swastika-like or otherwise symbol-like
  arrangements (blades with baffles, vortex patterns, crosses of jets);
- look at every rendered image before it is committed, at full size and as
  a thumbnail, and ask whether anyone could read a symbol into it; when in
  doubt, change the geometry, the time or the crop.

## README and documentation stay current

Every pull request that is merged to main leaves the README and the
documentation describing main as it is after the merge. When a pull
request adds, changes or removes a feature, it updates in the same pull
request:

- the full feature list [docs/features.md](docs/features.md) ("What works
  today"), and in the [README](README.md) the status line, the six-line
  summary "What works today" when a whole area changes, the commands under
  "Build and run" and the crate table when they change;
- the gallery entry, case files and recipe (above), the verification page
  of the stage, and `docs/keys.md` for new keys;
- the requirements (`docs/requirements.md`): the "Implemented" note of the
  decision or stage row;
- its changelog fragment.

A pull request that changes nothing a user can see says so in its
description. The checklist of the pull request template repeats these
points.

## Working in parallel

Several pull requests are often open at once. These rules keep them from
conflicting with each other, because GitHub runs no CI on a pull request
with a merge conflict, and every merge of main into a pull request runs
its CI again:

- **Changelog**: write the entry as a fragment `changelog.d/<PR number>-<short-name>.md`
  (see [changelog.d/README.md](changelog.d/README.md)), not in
  `CHANGELOG.md`; the release pull request collects the fragments.
- **Numbers**: gallery, recipe and case entries take the pull request's
  number (above); decisions in `docs/requirements.md` keep their D numbers,
  so add a decision in one pull request at a time.
- **Index pages**: regenerated after the merge (above).
- **Topic pages** of the gallery and recipes merge with git's `union`
  driver (`.gitattributes`): when two pull requests append to the same
  page, a local `git merge origin/main` keeps both entries without a
  conflict. Read the merged page once; GitHub itself ignores the driver
  and still reports the conflict until main is merged locally.
- **Do not merge main into a pull request whose CI is green**: merge the
  pull request. Merge main only for a real conflict.
- **The lock file of the Python crate** (`crates/foehncraft-py/Cargo.lock`)
  is separate, because the solver workspace builds without Python. After
  changing dependencies run
  `cargo update -w --manifest-path crates/foehncraft-py/Cargo.toml`; CI
  checks it first.
- **After a merge**, watch the main run of the merge commit (macOS and
  Windows are tested there, not in the pull request) and fix forward at
  once when it fails.
