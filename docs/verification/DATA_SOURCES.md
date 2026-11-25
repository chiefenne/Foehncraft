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
