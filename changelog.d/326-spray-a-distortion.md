### Added

- The Corrected Distortion model of Nguyen, Dahms, Pickett and Tagliante (2022) for distorted spray drops (ECN Spray A, SA4, D316): `[spray.breakup] drag = "corrected-distortion"`. A drop of TAB distortion y drags, heats and evaporates as an oblate or prolate spheroid, with the internal circulation of a viscous drop (Feng and Michaelides 2001) and the shape corrections of Richter and Nikrityuk (2012). A fully flattened drop drags 1.0 to 1.9 times a sphere instead of the 3.6 times of the dynamic drag, and evaporates 2 to 5 times faster. New case `spray_a_distortion.toml` (ECN Spray A, a quarter domain), gallery 326, recipe 326.

### Fixed

- With the Peng-Robinson surface, the saturation limiter of D253 gave gas above the critical mixing temperature (651 K at 6 MPa) no room for vapour, so it held back all the evaporation in the hot gas of a Spray A run. Drops then heated to 900 K without evaporating, and the liquid length grew with the vapour tip. Gas above the critical mixing temperature now takes vapour up to the critical surface fraction.
- `[spray.ecn]` liquid lengths of a half or quarter domain counted only the liquid inside the domain, so a line of sight across a symmetry plane saw half the spray. The projection now adds the mirror images of the parcels across the symmetry planes (and slip walls) that contain the injector axis. A sector whose images do not tile the whole spray is refused.
