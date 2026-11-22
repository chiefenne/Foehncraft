# Engineering Gallery

Industrial flows with complex geometry: stirred tanks, pumps, manifolds,
valves, ventilated rooms, heat sinks. Each entry answers an engineering
question (a power number, a pressure loss, a flow split) and shows how far
the answer can be trusted.

The [Example gallery](gallery.md) shows what the solver computes and
validates each feature against data. The Engineering Gallery uses those
features together on geometries where no reference data exist, so its
entries are trusted by checks that need no data (below), not by a
comparison with measurements. The proposal behind the track is
[step 25](requirements/step25-industrial-applications.md).

All commands run from the repository root. The cases are in
`examples/applications/<name>/`; each folder holds the mesh script, the
production case, a coarse smoke case for CI and the post-processing
script. The case files are listed in the
[Case catalogue](cases/applications.md), the generic workflows in the
[recipes](recipes/applications.md).

Each entry takes the number of the pull request that added it ("application
N"); its recipe and its case rows carry the same number. The numbers are
global and stay the same.

## How far to trust an entry

Every entry has a check card with the checks it ran and the checks it
skipped (a fine mesh may not fit the hardware of the run):

| Id | Check | Pass rule |
|---|---|---|
| C1 | Mass balance | net imbalance below 1e-6 of the throughput |
| C2 | Energy balance | boundary power = dissipation + numerical dissipation; numerical dissipation not negative |
| C3 | Iterative convergence | engineering quantities change by less than 0.1 % over the last 10 % of the iterations |
| C4 | Mesh sensitivity | three meshes with a ratio of at least 1.3; GCI below 5 % (or the band is reported) |
| C5 | Wall resolution | y⁺ in the range of the wall treatment on at least 95 % of the wall area |
| C6 | Mesh quality | no cell left unremediated; quality histograms in the entry |
| C7 | Plausibility | against an engineering correlation or handbook value where one exists; **not** validation |
| C8 | Independence of the domain | one extra run with a larger domain changes the quantity by less than 1 % |
| C9 | Health | no warning left open; every remaining warning explained |

Trust levels:

- **checked**: C1–C6 and C9 pass, C7 within the uncertainty of the
  correlation where one exists;
- **checked, model-sensitive**: the checks pass, but the result depends on
  the turbulence model by more than the GCI;
- **demonstration**: runs and looks right, checks incomplete (for example a
  single mesh); allowed for the first version of an entry.

## Pages

The track is split by industry sector; each page has its own table of
contents.

| Page | What it covers | Applications |
|---|---|---|
| [Rotating machinery](applications/rotating-machinery.md) | stirred tanks, pumps, fans | 131, 161 |
| [Heat exchangers](applications/heat-exchangers.md) | shell and tube, plate and finned exchangers (conjugate heat transfer) | 203 |

## All applications by number

| # | Application | Page |
|---|---|---|
| <a id="stirred-tank"></a>131 | [Stirred tank with a three-blade hydrofoil: power and pumping number](applications/rotating-machinery.md#stirred-tank) | [Rotating machinery](applications/rotating-machinery.md) |
| <a id="fyfe-stirred-tank"></a>161 | [Stirred 2 L lab tank of Fyfe et al.: power number and mixing time against measurements](applications/rotating-machinery.md#fyfe-stirred-tank) | [Rotating machinery](applications/rotating-machinery.md) |
| <a id="shell-tube-hx"></a>203 | [Shell and tube heat exchanger: a conformal multi-region tet mesh from a STEP assembly](applications/heat-exchangers.md#shell-tube-hx) | [Heat exchangers](applications/heat-exchangers.md) |
