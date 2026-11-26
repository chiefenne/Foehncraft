# Step 28: what ECN Spray A needs beyond breakup (D310 to D319)

Maintainer request of 2026-11-25: "yes, add all we need" for the Engine
Combustion Network (ECN) Spray A case, the non-reacting diesel reference
spray, and "check the papers in more depth, specifically the newest ones".
Companion of the [spray specification](step28-sprays.md) (D240 to D260) and
of the [validation page](step28-spray-validation-statistics.md) (D261 to
D268). The breakup and injection model (KH-RT with blob injection and a
rate-of-injection table) comes from the Spray G work and the named-liquid
property library (n-dodecane at low pressure) from the spray validation work;
this page covers everything else Spray A needs. Clean-room rules (D201)
apply: models from the open literature, measured data and our own
derivations only.

## 1. The case

ECN Spray A, inert variant (ECN5 Topic 3 guidelines, ECN Sandia pages, read
2026-11-25):

| Quantity | Value |
|---|---|
| Fuel | n-dodecane, 363 K at the nozzle |
| Nozzle | single axial hole; 210675: 0.0894 mm, Cd 0.90, Ca 0.98 (estimated), Cv 0.92; 210370: 0.0908 mm; 210677: 0.0837 mm |
| Injection | 150 MPa, 1.5 ms (short) or 5 to 5.5 ms (long); rate from the CMT injection-rate generator with the measured Cd, Ca |
| Ambient | 900 K, 22.8 kg/m³ (about 6 MPa), 0 % O₂, near quiescent (below 1 m/s) |
| Variants | 440, 700, 800, 1000, 1100 K; 15.2 and 7.6 kg/m³; 100 and 50 MPa |

Measured targets (Sandia, open download pages; licence terms not stated on
the pages, the published numbers are quoted with their citation and no data
set is stored, `DATA_SOURCES.md` rule 3):

- Liquid length: Mie 3 % threshold 10.0 mm (0.5 to 1.4 ms, bkldaAL1, 210677)
  and 10.4 mm (0.5 to 5.0 ms, bkldaAL4); projected liquid volume (PLV)
  threshold 0.2·10⁻³ mm³/mm², 11.2 mm quasi-steady (210370, 1 to 4 ms).
- Vapour penetration: schlieren, time series.
- Mixture fraction: pulse-burst Rayleigh, ensemble of 11 injections
  (Manin et al., ECN6), conversion Z = FA/(FA + 1), FA = (M_f/M_a) NFA.

Model definitions of ECN5 (appendix 2): liquid length is the farthest axial
position where the liquid volume fraction averaged over a 0.5 × 0.5 mm path
is 0.15 %; vapour penetration is the farthest axial position with fuel mass
fraction 0.1 %; fields are Favre averages on a uniform output grid (radial
0 to 20 mm by 0.04 mm, axial 0 to 100 mm by 0.2 mm); the smallest cell of
a submission is at most 0.25 mm.

## 2. What the newest literature says (read 2026-11-25)

- **The surface holds at 900 K.** Long-distance microscopy of n-dodecane
  droplets at Spray A conditions shows the classical regime (sharp surface,
  surface tension) at 897 K and 5.9 MPa (T_r 1.4, p_r 3.2), a transitional
  regime at 1056 K and 6.9 MPa and diffusive mixing at 1410 K and 9.2 MPa
  (Crua, Manin, Skeen, Pickett, SAE 2023-32-0117; Crua, Manin, Pickett, Fuel
  208, 2017). The transition is not instantaneous; the regime follows from
  the far-field pressure and temperature. A Lagrangian droplet with a sharp
  surface is therefore the right model for the 900 K target, but its surface
  state must be a high-pressure phase equilibrium, not Raoult's law.
- **Raoult's law fails at 6 MPa.** Peng-Robinson phase equilibrium at the
  surface raises the fuel vapour pressure and dissolves the ambient gas in
  the liquid (Luijten and Kurvers, Atomization and Sprays 20, 2010, abstract
  read: the effect is significant at diesel densities, the pressure effect on
  the mixture enthalpy is negligible). Our own result with the SA1 kernel
  (n-dodecane of Ambrose and Walton 1989 with the fitted α, N₂, k_ij =
  0.205 fitted to measured bubble points, VA.1 and VA.2): at 6 MPa the
  surface vapour mole fraction of the fuel is 1.5 to 1.9 times Raoult's
  value from 450 to 600 K, 8 to 20 mol % N₂ is dissolved in the surface
  liquid, the latent heat is 11 to 15 % below its value at saturation, and
  the two-phase surface ends at the critical mixing temperature 651.3 K
  (656.6 K at 3 MPa, 647.8 K at 8 MPa). He, Yi and Li (Appl. Therm. Eng. 245, 2024, abstract read) add a
  non-equilibrium term and report the equilibrium assumption as 8 % short on
  the spray liquid penetration above 3 MPa; noted, not built.
- **Equilibrium methods.** Peng-Robinson with van der Waals one-fluid mixing
  and a tangent-plane stability test (Matheis and Hickel, arXiv 1609.08533;
  Michelsen 1982), a TPn flash with Rachford-Rice (Fathi, Hickel, Roekaerts,
  arXiv 2205.07504, with the RKPR equation of state), real-fluid equilibrium
  at the droplet surface of a Lagrangian spray, X_s = X_l φ_l/φ_v with the
  rate law of Spalding unchanged (Srinivasan, Han, Poblador-Ibanez,
  Nocivelli, Yang, AIAA 2025-0742; volume-translated SRK, Chung et al. 1988
  transport, Takahashi's high-pressure diffusivity, the supercritical switch
  T > T_c). Perini et al. (SAE 2020-01-0240) solve the phase stability in a
  Lagrangian diesel spray and find the Raoult result off by up to 100 %.
  Collis et al. (arXiv 2604.23846, April 2026) model Spray A with a
  four-equation Eulerian method and attribute a slight over-prediction of the
  gas penetration to the ideal-gas law for the vapour.
- **Drop shape.** The Sandia "Corrected Distortion" model (Nguyen, Dahms,
  Pickett, Tagliante, Int. J. Multiphase Flow, 2022, doi
  10.1016/j.ijmultiphaseflow.2021.103927; built on Dahms and Oefelein, IJMF
  86, 2016) corrects drag and evaporation of distorted drops (viscous-sphere
  drag of Feng and Michaelides, distortion correction of Richter and
  Nikrityuk, sphericity correction of the Sherwood number, surface area of
  the spheroid). With the standard model the PLV liquid length is much too
  long; with the correction the liquid and mixture fields move towards the
  measurements.
- **Dense spray in a point-particle method.** Near the nozzle the liquid
  fills a large share of a 0.1 to 0.25 mm cell and the gas velocity a parcel
  sees is not resolved. Remedies in the literature: gas-jet theory for the
  relative velocity and a collision model (Abani et al., SAE 2008-01-0970,
  abstract read), exchange with a sphere of fixed physical radius instead of
  the cell (Torelli et al., Atomization and Sprays 2014), a volume-occupying
  dense core coupled to parcels (ELMO, Schmidt et al., arXiv 2107.03508);
  volume-filtered equations (Capecelatro and Desjardins 2013, already cited
  in step28-sprays.md 5.6). Collisions are often switched off for Spray A
  (Torelli et al.: little effect on the Sauter radius of an evaporating
  spray; the no-time-counter method makes mesh-aligned artefacts on
  Cartesian meshes, Schmidt and Rutland 2000).

Model results quoted in these papers (liquid lengths of their codes, their
mixing temperatures) are not used as targets (D201); only the measured ECN
numbers are.

## 3. Decisions

| ID | Decision |
|---|---|
| D310 | **Scope and staging** (section 5). Spray A gets, beyond breakup and injection: (1) real-fluid phase equilibrium and latent heat at the droplet surface, (2) high-pressure film properties, (3) a rule for the critical mixing state, (4) the ECN diagnostics, (5) dense-spray effects, (6) the case with its measured comparison. The breakup model, the blob injection and the rate table are not built here; the low-pressure liquid data (vapour pressure, density, c_l, latent heat, critical constants, vapour species data) come from the named-liquid library. |
| D311 | **Surface equilibrium: Peng-Robinson.** `[spray.evaporation] equilibrium = "raoult"` (default, today's model) or `"peng-robinson"`. Peng-Robinson (1976) with the 1978 κ(ω) for ω > 0.491 (Robinson and Peng 1978), van der Waals one-fluid mixing, a_ij = (1 − k_ij)(a_i a_j)^½, k_ij per pair as a key with a default fitted by us to measured equilibrium data (VA.2: 0.205 for N₂ + n-dodecane; the fuel's α by Mathias and Copeman 1983, fitted by us to the measured vapour pressure, VA.1; a k_ij of 0.30, the best value at 562 and 594 K alone, changes the surface vapour fraction at 6 MPa by 2 to 6 %). At the droplet temperature T_s and p₀, the fugacities of every species are equal in liquid and vapour; the inert gases keep their far-field ratios in the vapour. The surface vapour mass fraction Y_s of the fuel goes into the Spalding numbers unchanged. The dissolved gas shifts the equilibrium but is not booked: the liquid stays pure fuel in the mass and energy books (2 to 5 % by mass at the surface only, own estimate; recorded as a model limit). The solver is a fixed-count successive substitution from Wilson K values with a Newton finish, the same bits on CPU and GPU (D93 (6)); no convergence: refused for the sub-step with a counted diagnostic, Raoult used. |
| D312 | **Latent heat from the same equation of state.** With `"peng-robinson"` the heat to evaporate is the enthalpy of the fuel vapour as the gas books it, less the partial molar enthalpy of the fuel in the surface liquid from the Peng-Robinson departure function, so that gas and parcels close with one property set (D90 (3), D250). It tends to zero at the critical mixing state. Built as L(T, p₀) = L_lib(T) · L_PR(T, p₀)/L_PR,sat(T): L_PR = h̄_v − h̄_l of the fuel between the surface vapour and liquid (partial molar departure enthalpies), L_PR,sat the same for the pure fuel at its saturation pressure, L_lib the latent heat of the liquid data (Kirchhoff's law with the NASA-7 vapour). The gas receives the vapour at c_l T + L(T, p₀), the enthalpy the parcel gives up; with the ideal gas of today's gas solver this carries the departure enthalpy of the dense vapour (most of the 11 to 15 % at 6 MPa) into the gas as a sensible-heat deficit, so the books close exactly; a real-gas solver (D317) would carry it as departure enthalpy instead. |
| D313 | **High-pressure film properties.** Film mixture at the 1/3 reference state (Hubbard, Denny, Mills 1975): density from Peng-Robinson; diffusivity from the low-pressure value with a high-pressure correction (Takahashi 1975 or Riazi and Whitson 1993; chosen after both are read); viscosity and conductivity by Chung et al. (1988). Each source is read in full and gets its `DATA_SOURCES.md` entry before it is built. |
| D314 | **The critical mixing state.** When the surface equilibrium has no two-phase solution (T_s at or above the mixture critical temperature at p₀), the parcel is in the diffusive regime: no surface tension (no further breakup), no latent heat, and it mixes at the diffusion-limited rate of the d² law with the surface at the critical composition. Counted per step (parcels, mass) as a health number. At 900 K and 6 MPa this should not occur (section 2, own check); it matters for the 1000 and 1100 K variants. |
| D315 | **ECN diagnostics as tools.** Liquid length by the ECN5 definition (0.15 % path-averaged liquid volume fraction on 0.5 × 0.5 mm paths), by PLV (projected liquid volume per area from the parcels, threshold a key, default 0.2·10⁻³ mm³/mm²), and by liquid mass share (95, 97, 99 %); vapour penetration (fuel mass fraction 0.1 %); time- and Favre-averaged fields on the ECN output grid; axial and radial mixture-fraction profiles. Registered tools (D196), written to `result.json` and CSV with the definition in the header (D265). |
| D316 | **Dense-spray effects:** chosen by the maintainer (card, section 6). The card is open; the work goes on with its recommended option, distortion first, then the liquid volume in the gas equations and the spreading of the sources over a fixed radius (Torelli et al. 2014), collisions later. |
| D317 | **Real gas in the gas phase:** chosen by the maintainer (card, section 6). |
| D318 | **Turbulence for the first Spray A runs.** RANS k-ε with the existing case keys; the case states the value of `c_epsilon1` it uses and why (the round-jet spreading of standard k-ε is known to be too large; the value is set by a sensitivity run against the measured vapour penetration, not copied). Large-eddy simulation of Spray A needs 10⁷ cells and several realisations; later. |
| D319 | **Thermodynamics checked before the spray.** The equilibrium, the mixing line and a mixing-limited liquid-length model (Siebers 1999, an analytic scaling calibrated on measurements) are checked against measured data without CFD (VA.1 to VA.5), so a liquid-length error of the full case can be traced to mixing or to thermodynamics. |

## 4. Checks

| Check | Content | Passes when |
|---|---|---|
| VA.1 | Pure n-dodecane: Peng-Robinson saturation pressure with the fitted Mathias-Copeman α against the measured vapour-pressure equation of Ambrose and Walton (1989), valid to the critical point, 300 to 655 K | within 1 % (the plain κ(ω) of 1978 is 16 % off at 300 K) |
| VA.2 | N₂ + n-dodecane bubble pressures and vapour compositions against García-Córdova, Justo-García, García-Flores, García-Sánchez, J. Chem. Eng. Data 56 (2011) 1555, doi 10.1021/je1012372 (169 measured points, 344 to 594 K, up to 60 MPa; NIST ThermoML) | k_ij fitted by us; mean deviation of the bubble pressure quoted, at most 10 % up to 20 MPa |
| VA.3 | Low-pressure limit: the Peng-Robinson surface vapour fraction tends to Raoult's law as p → 0 (the Poynting factor and the second virial terms are physical and grow with p) | above Raoult by less than 0.5 % at 10 kPa and 3 % at 0.1 MPa, growing with p, 300 to 400 K |
| VA.4 | Adiabatic mixing of 363 K n-dodecane with 900 K N₂ at 6 MPa (own derivation): frozen and equilibrium mixing temperature over the mixture fraction, the saturation mixture fraction | the two-phase region and the end points exact to the energy balance (round-off). **Met (#320):** 7e-16 of the stream enthalpy difference on the frozen, equilibrium and ideal lines; Z_sat 0.357 (equilibrium) and 0.320 (ideal) |
| VA.5 | Siebers liquid length with our thermodynamics against the measured Spray A liquid lengths (900 K, and the temperature and density variants where the measured value is published) | ideal and real-fluid results quoted; the real-fluid one closer. **Not met (#320):** at 900, 700, 1200 K (22.8 kg/m³) and 900 K (7.6 kg/m³) the equilibrium Z_sat gives 7.2, 11.3, 4.9 and 10.7 mm, the ideal 8.2, 13.0, 5.7 and 11.3 mm, measured (DBI) 11.7, 16.3, 7.4 and 20.6 mm (C_L = 0.62, 21.5°, d = 90 µm, C_a = 0.95). The real fluid shortens the length by 5 to 14 %, away from the measurement; C_L was set with Siebers' own B, which is not open. The trends are right. Recorded as a gap of the SA1 feature row |
| VA.6 | Closed box with Peng-Robinson evaporation | mass and energy ledgers to round-off; pairs zero to the bit (D254) |
| VA.7 | Order independence with Peng-Robinson (V28.1) | bitwise on 1 and 4 threads |
| VA.8 | Spray A inert, 900 K, 22.8 kg/m³ | liquid length (ECN5 and PLV) against 10.4 and 11.2 mm, vapour penetration against schlieren, mixture fraction against Rayleigh; bands set with the case from the measurement uncertainty |
| VA.9 | Spray A temperature variants | the measured trend of liquid length with temperature |

## 5. Stages

| Stage | Content | Depends on |
|---|---|---|
| SA0 | This page | – |
| SA1 | Real-fluid kernel (Peng-Robinson, mixtures, fugacity, departure functions, surface equilibrium, stability test, mixing line), `equilibrium = "peng-robinson"`, D312 latent heat, D314 critical state; checks VA.1 to VA.4, VA.6, VA.7; gallery example with the equilibrium and mixing-line images. **Part 1 (#314):** kernel, surface table, `[spray.evaporation] equilibrium` and `[spray.evaporation.peng_robinson]` (`fuel`, or `critical_temperature`, `critical_pressure`, `acentric_factor`, `molar_mass`; `ambient`, `k_ij`, `table_step`), critical-state counts; VA.1 to VA.3, VA.6 (closed box, water in nitrogen); gallery 314. **Part 2 (#320):** Michelsen stability test, isothermal flash, the mixing line `realfluid::Streams` (frozen, equilibrium, ideal; Z_sat), the n-dodecane ideal gas (Scott 1974), Siebers' liquid length; VA.4, VA.5 (moved here from SA2: it needs only the mixing line); gallery 320. Still open: VA.7 with a fuel vapour species | the liquid library (n-dodecane vapour species) |
| SA2 | High-pressure film properties (D313); VA.5 moved to SA1 part 2. **Done (#324):** the real-fluid film of the Peng-Robinson surface (`[spray.evaporation.peng_robinson] film`, default `real-fluid`): Peng-Robinson density and c_p, Chung et al. (1988) dense-fluid factors on the kinetic-theory viscosity and conductivity, Riazi and Whitson (1993) on the diffusivity (Takahashi 1975 not built: a chart); NIST nitrogen within 1.3 % and 4.8 % on the factors; the surface table's early "critical" end from a cold start fixed; gallery 324. The coefficient tables are **to verify** against the printed sources | SA1, sources read |
| SA3 | ECN diagnostics (D315) | – |
| SA4 | Dense-spray effects (D316). **Part 1 (#326):** the Corrected Distortion model of Nguyen et al. (2022), `[spray.breakup] drag = "corrected-distortion"`: the drag, heat and mass transfer of the TAB spheroid with the internal circulation of a viscous drop (Feng and Michaelides 2001a, 2001b) and the shape ratios of Richter and Nikrityuk (2012), as factors on the sphere laws; gallery 326 with a first Spray A run (quarter domain, both drag laws). Also fixed there: the saturation limiter above the critical mixing temperature, and the `[spray.ecn]` projection of half and quarter domains (mirror images). Still open: the liquid volume in the gas equations and the source spreading (part 2); collisions | card (recommended option taken) |
| SA5 | Real gas in the gas phase if chosen (D317) | card |
| SA6 | The Spray A case: mesh, injection, runs at 900 K and the variants, VA.8, VA.9; gallery with liquid, vapour and mixture-fraction fields | breakup model and SA1 to SA4 |

## 6. Open points (maintainer cards)

- D316: which dense-spray effects to build first.
- D317: Peng-Robinson only at the droplet surface first, or also in the gas
  solver (density and enthalpy) now.

## 7. References of this page

- Abani, N. et al. (2008), An improved spray model for reducing numerical
  parameter dependencies in diesel engine CFD simulations, SAE 2008-01-0970
  (abstract read).
- Ambrose, D., Walton, J. (1989), Vapour pressures up to their critical
  temperatures of normal alkanes and 1-alkanols, Pure Appl. Chem. 61,
  1395-1403, doi 10.1351/pac198961081395 (read; eq. 3 and table 1 used
  for VA.1).
- Capecelatro, J., Desjardins, O. (2013): as in step28-sprays.md.
- Chung, T. H., Ajlan, M., Lee, L. L., Starling, K. E. (1988), Generalized
  multiparameter correlation for nonpolar and polar fluid transport
  properties, Ind. Eng. Chem. Res. 27, 671–679 (built in SA2 from the
  form printed in Poling et al. 2001 and secondary accounts, checked
  against NIST nitrogen; the printed paper **to verify**, DATA_SOURCES
  3g).
- Crua, C., Manin, J., Pickett, L. M. (2017), On the transcritical mixing of
  fuels at diesel engine conditions, Fuel 208, 535–548, doi
  10.1016/j.fuel.2017.06.091 (abstract and data page read).
- Crua, C., Manin, J., Skeen, S., Pickett, L. M. (2023), Transcritical mixing
  of fuels at reactive conditions, SAE 2023-32-0117 (read).
- Dahms, R. N., Oefelein, J. C. (2016), The significance of drop
  non-sphericity in sprays, Int. J. Multiphase Flow 86 (not yet read).
- ECN: ECN5 Topic 3 guidelines; Sandia Spray A data pages (liquid length
  bkldaAL1, bkldaAL4, nozzle 370 PLV; Rayleigh burst), read 2026-11-25.
- Feng, Z.-G., Michaelides, E. E. (2001a), Drag coefficients of viscous
  spheres at intermediate and high Reynolds numbers, J. Fluids Eng. 123,
  841-849; (2001b), Heat and mass transfer coefficients of viscous
  spheres, Int. J. Heat Mass Transfer 44, 4445-4454 (as printed by Nguyen
  et al. 2022; the originals **to verify**, DATA_SOURCES 3h).
- García-Córdova, T. et al. (2011), Vapor-liquid equilibrium data for the
  nitrogen + dodecane system at temperatures from (344 to 593) K and at
  pressures up to 60 MPa, J. Chem. Eng. Data 56, 1555–1564, doi
  10.1021/je1012372 (ThermoML record read).
- He, R., Yi, P., Li, T. (2024), A non-equilibrium evaporation model for a
  droplet and spray under diesel engine-relevant conditions, Appl. Therm.
  Eng. 245, 122855 (abstract read).
- Hubbard, G. L., Denny, V. E., Mills, A. F. (1975), Droplet evaporation:
  effects of transients and variable properties, Int. J. Heat Mass Transfer
  18, 1003–1008 (not yet read).
- Lemmon, E. W., Goodwin, A. R. H. (2000), Critical properties and vapor
  pressure equation for alkanes CnH2n+2, J. Phys. Chem. Ref. Data 29, 1-39,
  doi 10.1063/1.556054 (table 2 read: the normal boiling point check of
  VA.1).
- Liu, A. B., Mather, D., Reitz, R. D. (1993), Modeling the effects of
  drop drag and breakup on fuel sprays, SAE 930072 (the dynamic drag).
- Luijten, C. C. M., Kurvers, C. (2010), Real gas effects in mixing-limited
  diesel spray vaporization models, Atomization and Sprays 20, 595–609
  (abstract read).
- Matheis, J., Hickel, S., Multi-component vapor-liquid equilibrium model
  for LES and application to ECN Spray A, arXiv 1609.08533 (read).
- Mathias, P. M., Copeman, T. W. (1983), Extension of the Peng-Robinson
  equation of state to complex mixtures, Fluid Phase Equilib. 13, 91-108
  (the three-coefficient α; DOI to verify).
- Michelsen, M. L. (1982), The isothermal flash problem, Fluid Phase
  Equilib. 9 (not yet read).
- Nguyen, T. M., Dahms, R. N., Pickett, L. M., Tagliante, F. (2022), The
  Corrected Distortion model for Lagrangian spray simulation of
  transcritical fuel injection, Int. J. Multiphase Flow 148, 103927, doi
  10.1016/j.ijmultiphaseflow.2021.103927 (read in full; built in SA4).
- Peng, D.-Y., Robinson, D. B. (1976), A new two-constant equation of state,
  Ind. Eng. Chem. Fundam. 15, 59–64; Robinson, D. B., Peng, D.-Y. (1978),
  GPA research report RR-28 (the κ(ω) for heavy components).
- Perini, F. et al. (2020), An investigation of real-gas and multiphase
  effects on multicomponent diesel sprays, SAE 2020-01-0240 (read).
- Riazi, M. R., Whitson, C. H. (1993), Estimating diffusion coefficients of
  dense fluids, Ind. Eng. Chem. Res. 32, 3081–3088 (built in SA2 from
  secondary accounts; the printed paper **to verify**, DATA_SOURCES 3g).
- Richter, A., Nikrityuk, P. A. (2012), Drag forces and heat transfer
  coefficients for spherical, cuboidal and ellipsoidal particles in cross
  flow at sub-critical Reynolds numbers, Int. J. Heat Mass Transfer 55,
  1343-1354 (as printed by Nguyen et al. 2022; the original **to
  verify**).
- Schmidt, D. P., Rutland, C. J. (2000), A new droplet collision algorithm,
  J. Comput. Phys. 164 (cited by Torelli et al.).
- Schmidt, D. P. et al., The Eulerian Lagrangian Mixing-Oriented (ELMO)
  model, arXiv 2107.03508 (read).
- Siebers, D. L. (1999), Scaling liquid-phase fuel penetration in diesel
  sprays based on mixing-limited vaporization, SAE 1999-01-0528 (not yet
  read).
- Srinivasan, N., Han, A., Poblador-Ibanez, J., Nocivelli, L., Yang, S.
  (2025), Real-fluid equilibrium evaporation model for Lagrangian-Eulerian
  simulation of sustainable aviation fuel sprays at high-pressure
  gas-turbine conditions, AIAA 2025-0742 (read).
- Takahashi, S. (1975), Preparation of a generalized chart for the diffusion
  coefficients of gases at high pressures, J. Chem. Eng. Japan 7, 417–420
  (not built: the correlation is a chart; SA2 takes Riazi and Whitson).
- Torelli, R., D'Errico, G., Lucchini, T., Ikonomou, V., McDavid, R. M.
  (2014), A spherical volume interaction DDM approach for diesel spray
  modeling, Atomization and Sprays (read).
- Collis, H., Mirjalili, S., Khanwale, M., Mani, A., Iaccarino, G. (2026),
  An LES model with finite-rate phase change and subgrid spray based on a
  thermodynamically consistent four-equation multiphase model, arXiv
  2604.23846 (read in part).
