# Python

Part of the [Recipes](../recipes.md): the `foehncraft` Python package (D61): cases, expressions, runs with callbacks, results, materials and studies.

| # | Recipe |
|---|---|
| 33 | [Build, run and study cases from Python](#python-api) |

---

<a id="python-api"></a>
## 33. Build, run and study cases from Python

**Install.** The package builds its native module from the Rust sources
with maturin (Python 3.10 or newer, a Rust toolchain and CMake for HDF5).
Without it the package drives the `foehn` executable instead
(`backend="cli"`; `FOEHN=target/release/foehn` or `foehn` on `PATH`):

```
pip install ./python                     # or: cd python && maturin develop --release
pip install "./python[fields,study]"     # with h5py, NumPy and SciPy
python -c "import foehncraft; print(foehncraft.__version__, foehncraft.native_available())"
```

**Cases.** A `Case` is the table of a case file with dotted keys, the
names of `--set`; it writes TOML that reads back to the same table, and
keeps relative paths (mesh, tables) valid when saved elsewhere:

```python
from foehncraft import Case
from foehncraft.expr import y, t, constant, table, sin, pi

case = Case.load("examples/cases/flow/poiseuille_hex.toml")
case["fluid.viscosity"] = 0.05
case["constants"] = {"u_mean": 1.0}
case["boundary.xmin.velocity"] = [6 * constant("u_mean") * y * (1 - y), 0.0, 0.0]
del case["boundary.xmin.profile"]
case.save("output/my_case.toml")
```

**Expressions** are built with Python operators and render as the strings
of the case file (D53): `-x**2` is `-x ^ 2`, `(x < 1) & ~(y > 2)` is
`x < 1 and not y > 2`; `eq(a, b)`, `if_(c, a, b)`, `min_`, `max_`,
`clamp`, `step`, `table("ramp")(t)`. `foehncraft.expr.evaluate(e, x=…, y=…,
t=…, constants=…, tables=…)` evaluates with the solver's own compiler and
functions.

**Run** in process or through the command line; callbacks see every row
of the run record and a true return value stops the run (restart file,
final output, exit code 3):

```python
from foehncraft import run

res = run(case, out="output/my_case",
          sets={"solver.max_iterations": 2000},
          on_iteration=lambda r: r["continuity"] < 1e-8,   # steady runs
          on_step=None,                                     # transient runs
          log=print)
res.ok, res.status, res.exit_code, res.orders
res.value("errors.l2"), res.value("boundaries.ymin.force.0")
res.table("residuals")["continuity"]      # any CSV of the run as columns
res.fields()["velocity_x"]                # VTKHDF cell data (h5py)

from foehncraft.vtkhdf import internal     # the blocks of fields.vtkhdf in one piece
g = internal("output/my_case/fields.vtkhdf")
g["Points"], g["Connectivity"], g["Offsets"], g["CellData/pressure"]
```

**Materials** from the solver's database, and materials of a case:

```python
from foehncraft import materials
materials.get("water", temperature=313.15)["viscosity"]["value"]
case["materials.warm-oil"] = materials.material(
    base="engine-oil", viscosity=materials.power_law(0.25, 320.0, -6.0, source="own data"))
```

**Studies** (D62): variants from `grid`, `listing`, `latin_hypercube`
or `sobol` (the last two need SciPy), one case file and output directory
per variant, one table `study.csv`, resumed when run again:

```python
from foehncraft.study import Study, grid
rows = Study(case, grid({"fluid.viscosity": [0.05, 0.1, 0.2]}), "output/visc_study").run(jobs=1)
```

`jobs > 1` runs that many `foehn` processes at once. Examples:
`examples/python/poiseuille_inflow.py`,
`examples/studies/couette_viscosity_study.py` (gallery 24). Tests:
`python -m pytest python/tests`, lint `ruff check python`.
