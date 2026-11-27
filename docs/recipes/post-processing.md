# Post-processing

Part of the [Recipes](../recipes.md): derived fields with the output or afterwards from a restart file, and the statistics of transient runs. The keys are on the [key page](../cases/flow-keys.md#derived-fields-d72).

| # | Recipe |
|---|---|
| 30 | [Write vorticity, Q, λ₂, wall shear stress and other derived fields](#derived-fields) |
| 108 | [Write the Reynolds stress, the scales and the wall values of a turbulence model](#turbulence-derived-fields) |
| 110 | [Average a transient run: mean, RMS, covariance, confidence interval](#transient-statistics) |

---

<a id="derived-fields"></a>
## 30. Write vorticity, Q, λ₂, wall shear stress and other derived fields

**With the output.** List the fields in `[output] derived`; every
`fields*.vtkhdf` of the run then carries them (D72):

```toml
[output]
derived = ["vorticity", "q-criterion", "lambda2", "wall-shear-stress", "y-plus"]
```

The names: `vorticity`, `q-criterion`, `lambda2`, `enstrophy`,
`strain-rate`, `pressure-coefficient`, `total-pressure`,
`wall-shear-stress`, `y-plus`, and with a turning zone `relative-velocity`,
`frame-velocity`, `turbomachinery`, `relative-total-pressure`, `rothalpy`
([recipe 99](flow.md#rotor-fields)). The order of the list does not matter.
Cell fields go to the `internal` block; the wall shear stress and y⁺ go to
the boundary blocks (`boundary/<boundary>` in ParaView's block selector), zero on
boundaries that are not walls. The `units` attribute of the file lists every
dataset with its unit. The list may change during a run (re-read, D40).

**The pressure coefficient** needs a reference state; the density defaults
to `[fluid] density`, the velocity to that of `[output.coefficients]`:

```toml
[output.reference]
velocity = 1.0   # U_ref, m/s
pressure = 0.0   # p_ref, Pa
```

**After the run, from a restart file.** `foehn post derive` reads the
restart file and the case and writes the fields with the derived fields,
the file the run would have written at that iteration or step:

```
target/release/foehn flow examples/cases/flow/couette_derived.toml
target/release/foehn post derive examples/cases/flow/couette_derived.toml
# output/couette_derived/derived.vtkhdf (from restart.h5)

target/release/foehn post derive CASE.toml --restart output/CASE/restart_000400.h5 \
    --fields vorticity,q-criterion --file vortices.vtkhdf
target/release/foehn --json post derive CASE.toml   # summary as JSON
```

Without `--fields` it writes `[output] derived` of the case, or every
field the case allows when that list is empty (the pressure coefficient
only with a reference velocity, the turbulence fields with a model, the
wall heat values with `[energy]`, [recipe 346](energy.md#wall-heat)). `--restart` defaults to `latest`, the newest restart
file of the output directory; the file is `derived.vtkhdf` for
`restart.h5` and `derived_<N>.vtkhdf` for `restart_<N>.h5`. A run keeps
its restart file with `[run] restart_final = true` (or
`--set run.restart_final=true`).

**What they mean.** Q > 0 and λ₂ < 0 mark vortex cores (Hunt, Wray and
Moin 1988; Jeong and Hussain 1995); both are zero in pure shear, and Q is
also positive in boundary layers on convex walls. The wall values are those
of `wall_<boundary>.csv`: τ_w the tangential viscous stress of the fluid on
the wall, y⁺ = d √(τ_w/ρ)/ν with d the wall distance of the cell centroid.
Pictures: [gallery 22](../gallery/post-processing.md#derived-fields).

---

<a id="turbulence-derived-fields"></a>
## 108. Write the Reynolds stress, the scales and the wall values of a turbulence model

With a turbulence model (`[turbulence] model`), four more names go into
`[output] derived` or `foehn post derive --fields` (D72, specification
2.3):

```toml
[output]
derived = ["turbulence-scales", "reynolds-stress", "turbulence-wall"]
# the elliptic-blending k-ε also: "wall-normal"
```

- `turbulence-scales`: ε, ω, the time scale T, the length scale k^{3/2}/ε
  and ν_t/ν for every model. SST solves ω, so ε is β* k ω; the k-ε family
  solves ε, so ω is ε/(C_μ k). T is the one the model uses: 1/(β* ω) for
  SST, k/ε for the k-ε family, the bounded T of ζ-f and EB (also written
  as `turbulent_time_scale` with the fields).
- `reynolds-stress`: the six components of ⅔ k δᵢⱼ − 2 ν_t Sᵢⱼ, the
  stress the momentum equation sees (the eddy-viscosity models carry no
  other); compare `reynolds_stress_xy` with a measured −u'v'.
- `wall-normal`: n = ∇α/|∇α| of the elliptic-blending k-ε, the wall
  direction of its blending (0 where ∇α vanishes); refused with the other
  models.
- `turbulence-wall`: u_τ of the wall treatment and y* = d C_μ^¼ k_P^½/ν on
  the boundary blocks (zero off walls). y* is the wall coordinate of the
  k-ε wall functions and of the compound treatments; y⁺ is `y-plus`.
  `blending` is the viscous weight of the wall treatment: 1 where the wall
  is resolved (`low-re`, the sublayer of `standard`), 0 where the log law
  holds (the log layer, `scalable`), in between 1 − tanh((y⁺/5)²) of
  `automatic`, e^{−Γ} of `compound` and Γ = dU⁺/dy⁺ of `compound-2023`. A
  plot of it shows which walls of a mesh are resolved and which lean on
  the wall function. The
  heat-transfer coefficient and T⁺ of each scalar are the columns
  `<name>_h` and `<name>_t_plus` of `wall_<boundary>.csv`.

A laminar case refuses all four with CASE005. Pictures:
[gallery 108](../gallery/post-processing.md#turbulence-derived-fields).

---

<a id="transient-statistics"></a>
## 110. Average a transient run: mean, RMS, covariance, confidence interval

`[output.statistics]` (D71) averages a transient run in simulated time
from `start` to the end:

```toml
[output.statistics]
start = 5.0            # s; default the start of the run
batches = 10           # batches of the confidence interval (≥ 2, default 20)
phase_period = 0.3311  # optional: phase averages over this period (s)
phase_bins = 24        # bins per period (default 24)
fields = true          # default: also the fields (statistics.vtkhdf)
tolerance = 1e-3       # default: see "Is the window long enough?"
```

Every step is weighted with its time step, so an adaptive time step gives
time means, not step means.

- `statistics.vtkhdf` (with `[output] vtkhdf`): `mean_*` and `rms_*` of
  `velocity_x`, `_y`, `_z`, `pressure`, `k` (with a turbulence model) and
  every scalar by its name, the covariances `covariance_velocity_xy`,
  `_yz`, `_xz` (the resolved Reynolds shear stresses) and
  `covariance_velocity_x_<scalar>`, `_y_`, `_z_` (the resolved scalar
  fluxes). RMS is the standard deviation in
  time, √(⟨φ²⟩ − ⟨φ⟩²).
- `statistics.csv`: for every probe component `probes.<name>.u`, `.v`,
  `.w`, `.p` and the force coefficients `coefficients.cd`, `.cl`: mean,
  RMS, min, max, the half width of the 95 % confidence interval of the
  mean from `batches` batch means (Schmeiser 1982), the lag-1 correlation
  r1 of the batch means (von Neumann form) and a verdict.
- `statistics_phase.csv` (with `phase_period`): the time mean of every
  series in each bin of the phase (t − start)/P mod 1, the phase average
  of a periodic flow (a vortex street, a blade passing period).

**Is the window long enough?** The verdict is `sufficient` when the batch
means pass the von Neumann test of independence (Fishman 1978): their
lag-1 correlation r1 below 1.645 √((B − 2)/(B² − 1)), 0.47 for 10
batches, 0.35 for 20; else `too short`, and the run warns with STAT001. A
trend (a decay, a start-up still in the window) gives `too short`; start
later or run longer. With `phase_period` the batches hold whole periods
(as many as fit equally; the rest of the window enters the mean only), so
a periodic flow needs at least `batches` periods for the full number of
batches. A mean known well enough passes as well: when the half width,
widened for the correlation by √((1 + r1)/(1 − r1)), is below `tolerance`
(0.1 %) of the size √var + |mean| of the series, a slow drift far below
the size does not count (`tolerance = 0` keeps the test alone). The log
lists every series with its mean ± half width.

The section needs a transient run (CASE005) and checks `batches` and
`phase_period` (CASE004). Pictures:
[gallery 110](../gallery/post-processing.md#transient-statistics).
