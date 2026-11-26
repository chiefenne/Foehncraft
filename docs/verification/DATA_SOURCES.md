# Validation reference data: sources and attribution

The validation cases of Foehncraft compare against reference numbers from
published papers and published data sets. This page lists each source, the
quantities used, and how it is cited. **The repository holds only the numbers
that the tests and the validation pages quote, with their citation. It does
not hold the data sets themselves**, and it holds no figures, tables or prose
from the papers.

Status of each entry: **verified** means the citation (DOI or publisher page)
was found in a bibliographic search; **to verify** means the DOI, page range
or archive link still has to be checked against the paper or the archive
before the entry is used in a test.

The source code and the validation scripts are under the MIT licence (root
`LICENSE`). The reference values are factual scientific data points. If you
use them in your own work, cite the original authors and papers below, not
Foehncraft alone.

## 1. Plane channel flow, Re_τ = 395 and 590

### 1.1 DNS reference (Moser, Kim & Mansour 1999)

- Citation: R. D. Moser, J. Kim, N. N. Mansour, *Direct numerical simulation
  of turbulent channel flow up to Re_τ = 590*, Physics of Fluids 11(4),
  943–945 (1999).
- DOI: 10.1063/1.869966 (**verified**)
- Archive: the DNS data page of the authors' group,
  https://turbulence.oden.utexas.edu/MKM_1999.html (**verified**: the page
  names the friction Reynolds numbers 180, 395 and 590, cites Phys. Fluids
  11(4), page 943, and offers the statistics as a gzipped tar file and a
  browsable tree)
- Re_τ = 395 profiles: https://turbulence.oden.utexas.edu/data/MKM/chan395/profiles/
  (**verified**: the listing has chan395.flat, chan395.means,
  chan395.reystress, chan395.skew, chan395.velp and chan395.vortvar)
- Re_τ = 590 profiles: https://turbulence.oden.utexas.edu/data/MKM/chan590/profiles/
  (**verified**: the listing has chan590.flat, chan590.means, chan590.quad,
  chan590.reystress, chan590.skew, chan590.velp and chan590.vortvar). The 2015
  channel set (section 1.3) has no 395 or 590 case, so these values come from
  this archive
- Quantities used: mean streamwise velocity U⁺(y⁺); Reynolds stress
  components R₁₁⁺, R₂₂⁺, R₃₃⁺, R₁₂⁺; turbulent kinetic energy k⁺.

### 1.2 Model results of papers (not a data source)

The results of a model computed by its authors' own code (for example the
channel tables of an elliptic-blending model paper) are not reference
data (D201, rule 5 below). The numbering of 1.3 is kept.

### 1.3 Newer channel DNS set (Lee & Moser 2015)

- Citation: M. Lee, R. D. Moser, *Direct numerical simulation of turbulent
  channel flow up to Re_τ = 5200*, Journal of Fluid Mechanics 774, 395–415
  (2015).
- DOI: 10.1017/S0022112015002682 (**verified**: doi.org resolves to the
  publisher page for this article)
- Archive: the channel data page of the authors' group,
  https://turbulence.oden.utexas.edu/channel2015/ (the README asks that the
  data source be cited as above). The set has Re_τ = 180, 550, 1000, 2000 and
  5200. **It has no case at Re_τ = 395 or 590**, so the entries in section 1
  above keep their citation to Moser, Kim & Mansour (1999).
- Raw velocity fields are provided on request by the authors; statistics are
  published as text files on the archive. Nothing from the archive is copied
  into the repository. Only numbers quoted by a test or validation page, with
  this citation, may enter the repository.

## 2. Square duct (secondary flow, corner vortices)

### 2.1 Pinelli, Uhlmann, Sekimoto & Kawahara (2010)

- Citation: A. Pinelli, M. Uhlmann, A. Sekimoto, G. Kawahara, *Reynolds
  number dependence of mean flow structure in square duct turbulence*,
  Journal of Fluid Mechanics 644, 107–122 (2010).
- DOI: 10.1017/S0022112009992242 (**verified**: the publisher page carries
  this DOI, volume 644 and pages 107–122)
- Quantities used: mean streamwise velocity and the secondary velocity along
  the wall bisector and the corner diagonal, at the Reynolds numbers that
  the paper tabulates; the Re_b values to be taken from the paper (the
  abstract gives none).
- The data archive named in the draft (a KIT repository) has not been
  confirmed; it is removed until it is checked.

### 2.2 Gavrilakis (1992)

- Citation: S. Gavrilakis, *Numerical simulation of low-Reynolds-number
  turbulent flow through a straight square duct*, Journal of Fluid Mechanics
  244, 101–129 (1992).
- DOI: 10.1017/S0022112092002982 (**verified**: the publisher page carries
  this DOI, volume 244 and pages 101–129)
- Quantities used: secondary flow profiles and the wall shear distribution
  near the corners, as printed in the paper. The abstract states a bulk
  Reynolds number of 4410 (bulk velocity and hydraulic diameter); Re_τ is to
  be taken from the paper, not from this page.

## 3. Rotating channel

### 3.1 Kristoffersen & Andersson (1993)

- Citation: R. Kristoffersen, H. I. Andersson, *Direct simulations of
  low-Reynolds-number turbulent flow in a rotating channel*, Journal of Fluid
  Mechanics 256, 163–197 (1993).
- DOI: 10.1017/S0022112093002757 (**verified**, the publisher page carries
  this DOI)
- Quantities used: the asymmetric mean velocity profiles, the friction
  velocity ratio between the suction and the pressure side, and the
  Reynolds stresses across the channel, as printed in the paper for the
  rotation numbers it reports.
- The ERCOFTAC Classic Collection case is **not** used: its data carries
  CC BY-NC-SA 4.0, and it stays out of this repository. Values come from the
  paper.

## 3b. Stirred tank mixing (Fyfe et al. 2024)

- Citation: C. Fyfe, H. Barrington, C. M. Gordon, M. Reid, *A computer vision
  approach towards verifying CFD models*, Organic Process Research &
  Development 28 (9), 3661–3673 (2024), CC BY 4.0.
- DOI: 10.1021/acs.oprd.4c00229 (**verified**: the DOI and the page range
  are printed on the paper's pages)
- Quantities used: the measured plateau (macromixing) times of Table 3 for
  the phenolphthalein runs, 200 rpm unbaffled 17.3 s, 200 rpm baffled
  18.7 s, 400 rpm unbaffled 20.0 s, 400 rpm baffled 6.5 s (single runs,
  whole-reactor colour metric); the vessel, impeller and baffle dimensions
  of its Supporting Information. Not used: its CFD results (velocity planes,
  the Table 1 velocity differences) and the blend times of its
  correlation tool (Table 2), which are numbers of other solvers (D201).

### 3c. Power number and blend time correlations (to verify)

- Grenville blend time, θ₉₅ = 5.20 N_P^(−1/3) N⁻¹ (T/D)²: **to verify**
  (the source, Grenville's work on mixing times of turbulent stirred tanks,
  is quoted in the audit of 2026-11-21 as "Grenville (1992)"; no citation is
  confirmed yet).
- Measured power number of a four-blade 45° pitched-blade turbine in a
  baffled tank, about 1.2 to 1.4 in the turbulent range (Chapple, Kresta,
  Wall & Afacan 2002, Chem. Eng. Res. Des. 80): **to verify** before a test
  quotes it.

## 3d. Spray evaporation (SP4, SP5)

Published correlations, property formulations and analytic solutions; no
measured data set is stored. Each is quoted by the code or a test next to
its use.

- Saturation pressure of water: the auxiliary equation of the IAPWS-95
  formulation, Wagner, W., Pruss, A. (1993), International equations for
  the saturation properties of ordinary water substance, J. Phys. Chem. Ref.
  Data 22, 783-787, doi:10.1063/1.555926; the same equation as eq. (2.5) of
  Wagner, W., Pruss, A. (2002), The IAPWS formulation 1995 for the
  thermodynamic properties of ordinary water substance, J. Phys. Chem. Ref.
  Data 31, 387-535, doi:10.1063/1.1461829. The kernel test checks it at the
  triple point, 373.124 K and 101.325 kPa (values of the IAPWS-95 tables).
- Latent heat of water at 298.15 K, 2.4417 MJ/kg (IAPWS-95 tables, Wagner
  and Pruss 2002); at other temperatures from Kirchhoff's law with the
  vapour enthalpy of the NASA-7 data of GRI-Mech 3.0 (`H2O`, the species
  data the flow already uses) and h_l = c_l T.
- Droplet heat and mass transfer: Ranz, W. E., Marshall, W. R. (1952),
  Evaporation from drops, Chem. Eng. Prog. 48, 141-146 and 173-180
  (Nu = 2 + 0.6 Re^1/2 Pr^1/3, Sh = 2 + 0.6 Re^1/2 Sc^1/3); Abramzon, B.,
  Sirignano, W. A. (1989), Droplet vaporization model for spray combustion
  calculations, Int. J. Heat Mass Transfer 32, 1605-1618,
  doi:10.1016/0017-9310(89)90043-4 (film model with Spalding numbers B_M and
  B_T).
- d-squared law: Godsave, G. A. E. (1953), Studies of the combustion of
  drops in a fuel spray: the burning of single drops of fuel, Proc. 4th
  Symp. (Int.) on Combustion, 818-830, doi:10.1016/S0082-0784(53)80107-4;
  Spalding, D. B. (1953), The combustion of liquid fuels, Proc. 4th Symp.
  (Int.) on Combustion, 847-864, doi:10.1016/S0082-0784(53)80110-4. Used as
  the analytic solution of the single-droplet test (K = 8 ρ_g D ln(1 + B_M)
  / ρ_l at rest) and of the wet-bulb balance.
- Adiabatic saturation (closed box and V28.11): an exact energy and mass
  balance of the case itself, no external number.

## 3e. High-pressure droplet surface (Spray A, SA1)

Published measurements and equations, quoted by
`crates/foehncraft-particles/tests/high_pressure.rs` and the constants of
`foehncraft_kernels::realfluid`; no data set is stored.

- Vapour pressure of n-dodecane up to its critical point: eq. (3) and
  table 1 of Ambrose, D., Walton, J. (1989), Vapour pressures up to their
  critical temperatures of normal alkanes and 1-alkanols, Pure Appl. Chem.
  61, 1395-1403, doi:10.1351/pac198961081395 (T_c = 658.65 K, p_c =
  1.83 MPa, ω = 0.573, a = −9.08593, b = 2.77846, c = −5.1985,
  d = −4.173). The coefficients were read from the text layer of the open
  PDF; the test checks them against the normal boiling point 489.439 K of
  Lemmon, E. W., Goodwin, A. R. H. (2000), J. Phys. Chem. Ref. Data 29,
  1-39, doi:10.1063/1.556054, table 2 (0.1 % at 101.325 kPa). The
  Mathias-Copeman coefficients of `fuels::N_DODECANE` are our fit to this
  equation (300 to 655 K, every 5 K; reproduced by a test).
- Bubble points of nitrogen + n-dodecane: García-Córdova, T.,
  Justo-García, D. N., García-Flores, B. E., García-Sánchez, F. (2011),
  Vapor-liquid equilibrium data for the nitrogen + dodecane system at
  temperatures from (344 to 593) K and at pressures up to 60 MPa, J. Chem.
  Eng. Data 56, 1555-1564, doi:10.1021/je1012372 (record read through NIST
  ThermoML). The test quotes 21 of the points (three per isotherm, near 2,
  5 and 12 MPa); k_ij = 0.205 is our fit to all 75 points up to 20 MPa
  (mean deviation 2.7 %, at most 9.1 %).
- Critical points of the ambient gases from their reference equations of
  state (N₂: Span, Lemmon, Jacobsen, Wagner, Yokozeki 2000; O₂: Schmidt and
  Wagner 1985; CO₂: Span and Wagner 1996; H₂O: Wagner and Pruss 2002; Ar:
  Tegeler, Span, Wagner 1999); acentric factors from Poling, B. E.,
  Prausnitz, J. M., O'Connell, J. P. (2001), The Properties of Gases and
  Liquids, 5th ed., appendix A (**to verify** against the printed table).
- The Peng-Robinson constants Ω_a, Ω_b and Z_c are the roots of the
  critical conditions of the equation (Peng, D.-Y., Robinson, D. B. (1976),
  Ind. Eng. Chem. Fundam. 15, 59-64, doi:10.1021/i160057a011), recomputed
  to 19 digits by us; a test checks the triple root.

## 3f. Liquid property library and single droplets (D320, D321)

Only our fitted coefficients are in the repository
(`crates/foehncraft-kernels/src/liquids.rs`, `crates/foehncraft-flow/data/species-vapours.toml`);
the tables they are fitted to are not stored. Generators:
`tools/data/liquid_properties.py`, `tools/data/isooctane_table.py`,
`tools/data/vapour_species.py`.

- Water, n-heptane, n-dodecane: saturation tables of the NIST Chemistry
  WebBook, SRD 69, "Thermophysical Properties of Fluid Systems" (Lemmon,
  E. W., McLinden, M. O., Friend, D. G.), https://webbook.nist.gov/chemistry/fluid/,
  1 K steps from the triple point region to T_c (URLs printed by
  `liquid_properties.py --urls`). The tests quote single table points.
- Iso-octane (not in the NIST fluid system): vapour pressure, Wagner 3-6
  form of McGarry, J. (1983), Ind. Eng. Chem. Process Des. Dev. 22, 313
  (T_c 543.96 K, p_c 2.56155 MPa), checked against the Antoine law of
  Willingham, C. B., Taylor, W. J., Pignocco, J. M., Rossini, F. D. (1945),
  J. Res. NBS 35, 219 (297.5 to 373.3 K); density, specific heat, latent
  heat and viscosity: the DIPPR correlations printed in Perry's Chemical
  Engineers' Handbook, 8th ed. (2008), tables 2-32, 2-153, 2-150 and 2-313
  (published coefficients, transcribed through the `chemicals` Python
  package data files and **to verify** against the printed tables); latent
  heat checked against the calorimetric values of Majer, V., Svoboda, V.
  et al. (1979, 1982, 1985) quoted on the NIST WebBook (±0.1 kJ/mol);
  surface tension: Mulero, A., Cachadiña, I., Parra, M. I. (2012), J. Phys.
  Chem. Ref. Data 41, 043105 (σ0 = 0.04794 N/m, n0 = 1.209, 273 to 373 K);
  normal boiling point 372.4 K (NIST WebBook).
- Vapour species NC7H16 and IC8H18: ideal-gas heat capacities of the NIST
  WebBook gas-phase data (n-heptane and 2,2,4-trimethylpentane pages) and
  their enthalpies of formation (−187.8 and −224.1 kJ/mol); Lennard-Jones
  parameters from the critical point and acentric factor by Tee, L. S.,
  Gotoh, S., Stewart, W. E. (1966), Ind. Eng. Chem. Fundam. 5, 356-363.
- Single droplets: Chauveau, C., Birouk, M., Halter, F., Gökalp, I. (2019),
  Int. J. Heat Mass Transfer (title, volume and pages **to verify**; the
  rates are those of its eq. (4) as read in this work) (microgravity
  K0 = 3.6552e-4 T − 0.1078 mm²/s, 473 to 748 K, n-heptane in N2 at
  0.1 MPa, 14 µm cross fibres, d0 about 0.37 mm); Chauveau, C., Halter, F.,
  Lalonde, A., Gökalp, I. (2008), ILASS-Europe, Como (fibre-free
  extrapolation K0 = 4.485e-4 T − 0.1369 mm²/s). Measurements, D201
  accepted class. The model results of Both, A., Mira, D., Lehmkuhl, O.
  (2022), Int. J. Heat Mass Transfer, arXiv:2201.02818, are quoted only as
  context for the gap, not as targets (D201).
- Film reference state: Hubbard, G. L., Denny, V. E., Mills, A. F. (1975),
  Int. J. Heat Mass Transfer 18, 1003-1008; Yuen, M. C., Chen, L. W. (1976),
  Combust. Sci. Technol. 14, 147-154.

## 3g. ECN Spray G, condition G1 (D309)

Measurements of the Engine Combustion Network, quoted by
`tools/plot/spray_g_g1.py` (the penetrations) and
`examples/cases/flow/spray_g_g1_injectors.toml` (the injector); the files
are not stored.

- Liquid and vapour penetration along the injector axis, iso-octane, G1
  (573 K, 3.5 kg/m³, 20 MPa): the Sandia Spray G data of the ECN,
  https://ecn.sandia.gov/data/sandia-spray-g-data/ (files
  `liquid_data_2e-4_thr.mat`, `liquid_data_2e-3_thr.mat`, column 16, and
  `vapor_data.mat`, column 12; high-speed extinction imaging evaluated as
  projected liquid volume above 0.2e-3 and 2e-3 mm³/mm², and schlieren
  imaging, by J. Hwang and L. M. Pickett), cited as the data page asks:
  Hwang, J., Weiss, L., Karathanassis, I. K., Koukouvinis, P., Pickett,
  L. M., Skeen, S. A. (2020), Spatio-temporal identification of plume
  dynamics by 3D computed tomography using engine combustion network
  spray G injector and various fuels, Fuel 280, 118359,
  doi:10.1016/j.fuel.2020.118359. The script quotes 12, 10 and 12 points
  (the frames nearest to 0.05 to 1.2 ms). Measurements, D201 accepted
  class.
- Plume direction 33°: our mean of `plume_angles.mat` of the G1 3D
  computed tomography of the same data page and paper (0.3 to 0.7 ms,
  planes 2 to 15 mm from the tip: 32.9° to 33.7°).
- Rate of injection: the ECN file `roi_20MPa_680us.csv` (20 MPa, 680 µs
  energising, 10.16 mg); the injector file holds our 17-point
  least-squares fit to it, not the file.
- Nozzle: ECN Spray G specification (8 holes of 0.165 mm, orifice length
  0.16 to 0.18 mm, drill angle 37°), the ECN page "Spray G geometry".
  Area contraction 0.635 to 0.687 measured by Payri, R., Gimeno, J.,
  Martí-Aldaraví, P., Vaquerizo, D. (2015), Momentum flux measurements on
  an ECN GDi injector, SAE 2015-01-1893, doi:10.4271/2015-01-1893, and the
  full plume cone of 25° from the tomography and the schlieren boundary
  at 50 % intermittency, both as quoted by Colmenares F., J. D., Ameen,
  M. M., Patel, S. S., Large eddy simulation of gasoline sprays in a
  Lagrangian-Eulerian framework using the high-order spectral element
  method, ASME J. Eng. Gas Turbines Power, paper GTP-22-1576 (OSTI
  2427994; volume and pages **to verify**). Their model settings (B₁ = 80
  and others) are not used (D201).
- Drop size: the phase-Doppler data of the same data page ("Droplet size
  and velocity at z = 15 mm", GM, S. Parrish, injector #16, plume 1, file
  `PDI_distribution.xlsx`): the time-averaged Sauter mean diameter of
  9 µm, quoted by the plot script. Measurements, D201 accepted class.
- Turbulence scale of KH-ACT at low gas density (a finding quoted in the
  gallery text, no numbers taken): Magnotti, G. M. (2017), Modeling the
  influence of nozzle-generated turbulence on diesel sprays, PhD thesis,
  Georgia Institute of Technology (OSTI 1785739), section 1.2.3.3.

## 4. Rules for adding an entry

1. A value goes into the repository only when a test or a validation page
   quotes it, with its citation next to it.
2. Digitised values from a figure are marked as such and name the figure.
3. Full data sets, archives and fetch scripts are not added to the
   repository. Data needed for a validation run is downloaded in the
   project for that run only, and stays out of the repository.
4. A new entry needs a verified DOI or publisher page, or it is marked
   **to verify** until it is.
5. Reference numbers are measurements, direct numerical simulations
   (published scalar values with their citation, no data sets), published
   correlations and analytic solutions. Numbers computed by another
   RANS or CFD code, whether a solver's result on a benchmark page or the
   authors' own model result in a paper, are not used as targets (D201,
   decision of 2026-11-22). Direct numerical simulations stay an accepted
   class: D201 targets code-to-code RANS/CFD numbers, not DNS.
