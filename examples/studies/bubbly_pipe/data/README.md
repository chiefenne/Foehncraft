# Bubbly up-flow in a vertical pipe: measured void profiles

`wang1987_void.csv` holds the measured radial profiles of the void
fraction of the bubbly pipe flow of Wang et al. (1987): air and water in
a vertical pipe of D = 57.15 mm, up-flow with the liquid superficial
velocity j_L = 1.08 m/s and the gas superficial velocities
j_G = 0.1 and 0.4 m/s, measured with hot-film probes at L/D = 35. Both
profiles peak near the wall (wall-peaked bubbly flow).

| Column | Content |
|---|---|
| `jg` | gas superficial velocity j_G (m/s) |
| `r_over_R` | radial position r/R |
| `alpha` | measured void fraction |

**Source.** The original paper is not open; the points were digitised
(2026-11-28) from Fig. 2 (left panel, "measurement" symbols) of
W. K. In and T. H. Chun, "CFD Simulations of Two-Phase Bubbly Flows in a
Circular Channel and a Rod Bundle Subchannel", Korea Atomic Energy
Research Institute, INIS/ETDEWEB record FR0200356,
https://www.osti.gov/etdeweb/servlets/purl/20241563, which reproduces the
measurements of Wang et al. The reading uncertainty is about ±0.005 in
α and ±0.01 in r/R; the steep wall peak (r/R > 0.9) has fewer points
than the original.

**Reference.** S. K. Wang, S. J. Lee, O. C. Jones Jr. and
R. T. Lahey Jr., "3-D turbulence structure and phase distribution
measurements in bubbly two-phase flows", International Journal of
Multiphase Flow 13 (3), 1987, pp. 327–343. The bubble diameter of
3 mm of the simulations is not stated in the original; it follows
Colombo & Fairweather, International Journal of Multiphase Flow 77
(2015), p. 222 ff.

**Licence.** Measured values are facts; the digitised table is an
independent reading of a published figure and is redistributed with
this repository (MIT). Cite Wang et al. (1987) when you use it.

Used by `examples/studies/bubbly_pipe/compare.py` (stage MP3, D356).
