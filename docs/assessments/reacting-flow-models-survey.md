# Reacting flow models after ECFM-3Z: survey and plan (D360 to D365)

Status: survey 2026-11-28, decisions proposed, **no model coding before the
maintainer's go**. Requested by the maintainer (2026-11-28): "what is the latest
and greatest compared to the ECFM-3Z model", then "research and card me".

Clean-room (D201): only published papers, public workshop summaries and the
public manuals of codes are used; every number below carries its source.
Commercial codes are named where they show who uses a model (rule of
2026-11-22). Numbers marked *(figure/summary)* were read from a figure or a
machine summary of a long paper and should be re-read before they become a
pass band.

## 0. Short answer

ECFM-3Z (Colin and Benkenida, 2004) is a flame-surface model for engines: it
tracks how much flame surface there is and how fast it burns, with three mixing
zones (pure fuel, pure air with residual gas, mixed) so that it can do diesel
diffusion burning as well as spark-ignition fronts. It is still a solid RANS
engine model, and its own family moved on to ECFM-CLEH (four zones, with a
chemical-equilibrium limit of the burn rate, Abouri et al., 2008) and ECFM-LES
(Richard et al., 2007).

What leading groups use today instead, ranked by how the evidence reads:

1. **Tabulated flamelets** (FGM, flamelet/progress variable FPV, unsteady FPV,
   the tabulated flamelet model TFM). Best accuracy per core hour on ECN Spray A
   (ignition delay within about 5 to 20 %, lift-off within about 4 to 14 %, at
   1/40 of the cost of detailed chemistry per cell) and the standard for gas
   turbines and premixed LES. Weakness: extinction and re-ignition, and lean
   hydrogen unless the table carries stretch and preferential diffusion.
2. **Detailed chemistry per cell with acceleration** (called "well-mixed" or
   "no turbulence–chemistry model": CONVERGE SAGE with multi-zone, Forte with
   dynamic adaptive chemistry and cell clustering, OpenFOAM with TDAC/ISAT).
   The default of engine codes because it needs no flame-type choice, but on
   Spray A it gives ignition delay and lift-off that are too long (up to 67 %)
   because it ignores the turbulent fluctuations in each cell.
3. **Finite-rate with a reactor closure** (EDC, PaSR). The furnace workhorse;
   the standard EDC constants fail in MILD and low-Damköhler flames, fixed by
   local constants (Parente et al., 2016) or PaSR.
4. **Transported PDF** (Lagrangian Monte Carlo or Eulerian stochastic fields)
   and **CMC / MMC**. The most accurate on lift-off in the ECN comparisons, and
   the most expensive (tens of times the cost of tables).
5. **Premixed-front models** for spark ignition (G-equation, thickened flame,
   flame surface density as in ECFM-LES). Still used in SI engine LES; none is a
   clear winner.
6. **Machine learning** (neural-network tables, learned chemistry): research
   only in 2026; no production use documented.

Our position: SF2b already does detailed chemistry per cell (laminar, CPU and
GPU), SF3 builds EDC now, SF13 tables are planned. The recommended next model is
the **non-adiabatic flamelet table (FGM/FPV, SF13 brought forward)** with one
table kernel for the furnace and for reacting sprays (card 1, D360).

## 1. ECFM-3Z and its family

| Model | Source | What it is |
|---|---|---|
| ECFM | Colin et al., OGST 58 (2003) | Coherent flame model with a transported flame surface density, for premixed and partially premixed SI engines |
| **ECFM-3Z** | O. Colin, A. Benkenida, "The 3-Zones Extended Coherent Flame Model (ECFM3Z) for computing premixed/diffusion combustion", Oil & Gas Sci. Technol. 59(6) (2004) 593–609, doi 10.2516/ogst:2004043 ([abstract](https://ogst.ifpenergiesnouvelles.fr/articles/ogst/abs/2004/06/colin_vol59n6/colin_vol59n6.html)) | ECFM plus three mixing zones per cell: pure fuel, pure air with possible residual gas, and a mixed zone in which the ECFM flame model acts (abstract). Built to extend ECFM to unmixed diesel burning; auto-ignition modelled separately. Heat release is reported per mode: auto-ignition, flame propagation, post-flame oxidation, spark ignition ([STAR-CCM+ reference](https://stevedocs.azurewebsites.net/GUID-73CD350D-9A97-4BEB-B252-347894F4487C.html)) |
| **ECFM-CLEH** | Abouri, Duranti, Ravet, Zellat (2008), "Advances in combustion modeling in STAR-CD: validation of ECFM CLE-H model to engine analysis" ([ERC record](https://erc.wisc.edu/publications/Advances-in-Combustion-Modeling-in-Star-Cd-Validation-of-Ecfm-Cle-H-Model-to-Engine-analysis/)) | Four zones (premixed, diffusion, post-oxidation, auto-ignition inside the premixed zone), a premixed-to-diffusion transfer and burn-rate limits towards chemical equilibrium ([STAR-CCM+ reference](https://stevedocs.azurewebsites.net/GUID-3CFE0872-ADE8-42E5-A19F-249EC3BB99EA.html)) |
| **ECFM-LES** | Richard, Colin, Vermorel, Benkenida, Angelberger, Veynante, Proc. Combust. Inst. 31 (2007) 3059–3066, doi 10.1016/j.proci.2006.07.086 | Flame surface density transport for LES of SI engines; used for cycle-to-cycle variation with the ISSIM-LES spark model (D'Adamo et al., SAE Int. J. Engines 12(1), 2019, [pdf](https://iris.unimore.it/bitstream/11380/1171574/4/REV%20R1-Understanding%20the%20Origin%20of%20CCV%20using%20LES.pdf)) |

Who ships it: CONVERGE ([ECFM3Z from IFPEN](https://convergecfd.com/ecfm3z)),
STAR-CD/STAR-CCM+ (ECFM-3Z, ECFM-CLEH), Ansys CFX (ECFM).

Flaws we must not copy: the flame surface model needs a laminar flame speed
correlation and a wrinkling closure per fuel, a separate auto-ignition model and
a separate soot and NOx model; each has its own constants that are tuned per
engine. The public manuals list five burn-rate limit coefficients for CLEH
(BtPremLim, Ktrans, KtCleh, GapCmin, KtPstx) without the formulas. For a new
fuel (hydrogen, ammonia, e-fuels) every one of those has to be found again.
That is the main reason the field moved to chemistry-based models.

## 2. The model families today

### 2.1 Detailed chemistry per cell, with acceleration

Every cell is a small reactor integrated with the full mechanism over the time
step; the turbulence only mixes. Cheap to set up, works for any fuel and
regime, and the reason the main engine codes made it their default.

Evidence:
- On Spray A it ignores the sub-grid fluctuations: well-stirred reactor URANS
  gives "consistently longer ignition delays and lift-off lengths", up to 67 %
  apart from the stochastic-field PDF result, which stays within 14 % of
  experiment (Pang et al., Combust. Flame, 2018,
  [pdf](https://pure-oai.bham.ac.uk/ws/portalfiles/portal/48293044/CNF_D_17_00146R2_Manus_Figs.pdf)) *(summary)*.
  ECN2 (2012): well-mixed models "over-predicted ignition delay and lift-off
  length" ([ECN2 summary](https://ecn.sandia.gov/workshop/ECN2/ECN2proceedingsSummary.pdf)) *(summary)*.
- Cost: on Spray A/D in CONVERGE, well-mixed SAGE about 50.8·10³ core-hours
  against 1.2·10³ for unsteady FPV, about 40 times more (García-Oliver et al.,
  2020, [ScienceDirect](https://www.sciencedirect.com/science/article/pii/S2666691X20300282)) *(summary)*.
- Maghbouli et al., SAE 2016-01-0577 (heavy-duty optical engine, Spray B):
  "the flamelet model better predicts the sensitivity of lift-off length while
  the well-mixed model better predicts ignition delay"
  ([SAE](https://saemobilus.sae.org/content/2016-01-0577)).

Acceleration (what makes it usable):
- **Multi-zone / cell clustering**: cells grouped by temperature and
  equivalence ratio, chemistry solved per group and mapped back. CONVERGE
  licensed the LLNL multi-zone model: "same accuracy while reducing the
  simulation time by a factor of 10" (Flowers, DEER 2012,
  [slides](https://www.energy.gov/sites/prod/files/2014/03/f8/deer12_flowers.pdf);
  Raju et al., SAE 2012-01-0135).
- **Dynamic adaptive chemistry**: per-cell or per-cluster DRGEP reduction each
  step (Ansys Forte, [manual 5.3.2](https://ansyshelp.ansys.com/public/Views/Secured/corp/v252/en/forte_th/i85722.html)).
- **Cell agglomeration + pre-reduced mechanisms (SPARC)**: up to 41.5× on the
  chemistry and 13.3× overall on a sooting ethylene flame with a 185-species
  mechanism, temperature error 1·10⁻⁴, soot 1.5 % (Cuoci, Nobili, Parente,
  Grenga, Pitsch, Proc. Combust. Inst. 40 (2024) 105386,
  [pdf](https://re.public.polimi.it/bitstream/11311/1278933/1/1-s2.0-S1540748924001949-PROCI24.pdf)).
- ISAT (Pope 1997) and TDAC (Contino et al., PROCI 33, 2011) in OpenFOAM.
- GPU and better linear algebra: the LLNL adaptive-preconditioned solver cut a
  857-species chemistry time about 26× (DEER 2012 above).

### 2.2 Finite-rate chemistry with a reactor closure: EDC and PaSR

The cell is split into fine structures that react and surroundings that do not
(EDC, Magnussen 1981, 2005), or the rate is scaled by the ratio of chemical to
mixing time (PaSR). Used in every general-purpose code for furnaces, and what
SF3 builds now.

Evidence:
- Standard constants (C_γ = 2.1377, C_τ = 0.4083) over-predict the peak
  temperature in MILD / jet-in-hot-coflow flames and are said not to apply
  below Re_t ≈ 64 (Parente, Malik, Contino, Cuoci, Dally, Fuel 163 (2016)
  98–111, [pdf](https://re.public.polimi.it/bitstream/11311/1003924/4/Parente_et_al_2016.pdf)).
  Adelaide JHC HM1 at z = 120 mm: T_max error 27.8 % (standard) → 9.35 % with
  local C_τ ∝ 1/(Da_η √(Re_t+1)), C_γ ∝ (Da_η (Re_t+1))^½. Fixed changes by De
  et al. (C_τ = 3) and Evans et al. (C_τ = 3, C_γ = 1) do part of this.
- Bösenhofer, Wartha, Jordan, Harasek, Energies 11 (2018) 1902 compare the
  fine-structure treatments EDC 1981, 1996, 2005, 2016 and 2017; they agree at
  standard constants in strong turbulence
  ([pdf](https://repositum.tuwien.at/bitstream/20.500.12708/20063/1/Boesenhofer-2018-Energies-vor.pdf)).
- 20 kW quasi-industrial MILD furnace: EDC "fails in providing a reasonable
  estimation of the ignition region"; global temperature error EDC-2016 7.4 %,
  PaSR (C_mix = 0.7) 4.5 %, PaSR with a dynamic mixing time 2.3 %; outlet O₂
  3.32 against 3.30 % measured (Ferrarotti et al., Energy Fuels 32 (2018)
  10228, [pdf](https://orbi.umons.ac.be/bitstream/20.500.12907/37406/1/Key%20Modeling%20Aspects%20in%20the%20Simulation%20of%20a%20Quasi-industrial%2020%20kW.pdf)) *(part summary)*.
- LES of the Adelaide JHC: PaSR centreline peak +15 K, laminar finite rate
  −40 K (Li, Cuoci, Parente, Proc. Combust. Inst. 37 (2019) 4519,
  [pdf](https://re.public.polimi.it/retrieve/e0c31c0e-e619-4599-e053-1705fe0aef77/1-s2.0-S1540748918306382-main_Li_PROCI.pdf)).
- Sandia Flame D, RANS: EDC centreline peak temperature +500 K at x/D = 15 and
  NO about three times too high; an unsteady flamelet post-step gives +150 K
  and twice (Saini, Prakash, De, Yadav, [arXiv 2103.01450](https://arxiv.org/pdf/2103.01450)).
- Ammonia: RANS-EDC and LES-PaSR are among the common choices; EDC
  over-predicts NO in MILD furnaces (Vigueras-Zuniga et al., Renew. Sustain.
  Energy Rev. 244 (2027) 117566, [pdf](https://orca.cardiff.ac.uk/id/eprint/190194/1/1-s2.0-S1364032126008683-main.pdf)).

### 2.3 Tabulated flamelets: FGM, FPV, UFPV, TFM

Laminar flames (counterflow diffusion flames, premixed flames, or igniting
unsteady flamelets) are computed once with the full mechanism; the flow carries
a few control variables (mixture fraction, progress variable, their variances,
enthalpy deficit, scalar dissipation) and looks the rest up. Fluctuations enter
through presumed PDFs. Fuel-flexible as long as the table is rebuilt.

Evidence on ECN Spray A (measured ID 0.40–0.44 ms, lift-off 16.5–17.7 mm across
four labs, [ECN2](https://ecn.sandia.gov/workshop/ECN2/ECN2proceedingsSummary.pdf)):
- RANS FGM (Lib-ICE/OpenFOAM, reduced Yao mechanism): ID 0.91/0.45/0.25 ms
  against 0.85/0.41/0.24 ms at 800/900/1000 K; lift-off 25.1/16.5/11.7 mm
  against 26.2/16.7/12.2 mm (Maghbouli et al., Combust. Theory Model. 2018,
  [pdf](https://re.public.polimi.it/retrieve/e0c31c0c-d46a-4599-e053-1705fe0aef77/MaghbouliCTM2018.pdf)) *(summary)*,
  so about 5 to 10 % on ID and 1 to 4 % on lift-off.
- LES FGM (Aalto, 130 species, 11.5 M cells): ID 10–20 % short, lift-off 4–14 %
  off at 50/100/150 MPa (Kahila et al., Combust. Flame 191 (2018) 142,
  [pdf](https://acris.aalto.fi/ws/portalfiles/portal/17669459/1_s2.0_S0010218018300087_main.pdf)).
- RANS unsteady FPV in CONVERGE: ID and lift-off "in excellent agreement" at
  800–900 K, 40× cheaper than well-mixed (García-Oliver et al. 2020, above).
- Non-adiabatic 5D FGM: baseline ID 0.4525 against 0.435 ms (Zhang et al.,
  Combust. Theory Model. 2019, [pdf](https://core.ac.uk/download/326218187.pdf)) *(summary)*.
- The TFM of Argonne/CONVERGE (Kundu et al.) uses flamelet residence time
  instead of a progress variable, with a 2755-species mechanism (cited in the
  two papers above). The primary OSTI paper did not load; not verified here.

Weaknesses:
- Extinction and re-ignition: FPV is better than finite-rate LES on Sandia
  Flame D but "significantly worse" on Flame F because it over-predicts
  extinction (Yang et al., Combust. Sci. Technol. 2020,
  [pdf](https://www.yang.gatech.edu/publications/Journal/CST%20(2020,%20Sandia%20flames).pdf)).
- Lean hydrogen: unstretched premixed manifolds over-predict the consumption
  speed by about 20 %; a fixed-strain flamelet stays under 5 % (Porcarelli et
  al., 2025, [arXiv 2510.03324](https://arxiv.org/html/2510.03324v1)). In
  thermodiffusively unstable H₂ flames the H radical must be in the enthalpy
  cross-diffusion term of the table (Schepers and van Oijen, 2025,
  [arXiv 2503.02620](https://arxiv.org/pdf/2503.02620)).
- Ammonia: adiabatic FGM fails on NO and CO; non-adiabatic FGM is needed
  (Honzawa et al., in the 2027 review above). LES-FGM with extra NO and NO₂
  transport equations predicts emissions within 10 % (RANS 35 % low) on a
  swirl NH₃/H₂ burner (Mazzotta et al., Flow Turb. Combust. 114 (2025) 1043,
  [Springer](https://link.springer.com/article/10.1007/s10494-025-00638-7)).
- Memory: 5D tables reach 1.5 to 17.5 GB; tensor-train compression cuts them
  100 to 800× with an a-priori error bound and lookups up to 2.4× faster
  (Demir et al., 2026, [arXiv 2603.20240](https://arxiv.org/pdf/2603.20240)).
  Good for GPU memory.
- Review: van Oijen et al., "State-of-the-art in premixed combustion modeling
  using flamelet generated manifolds", Prog. Energy Combust. Sci. 57 (2016) 30.

### 2.4 Transported PDF and CMC/MMC

The joint PDF of composition is transported (Lagrangian particles or Eulerian
stochastic fields) so chemistry is closed exactly and only mixing is modelled;
CMC solves conditional equations in mixture fraction space.

- ECN2: CMC (ETH Zürich) "greatly improved" lift-off and slightly improved
  ignition delay over well-mixed; transported PDF (UNSW, Penn State) similar;
  cost "significantly larger" than well-mixed; flamelets "the most efficient"
  ([ECN2](https://ecn.sandia.gov/workshop/ECN2/ECN2proceedingsSummary.pdf)) *(summary)*.
- Eulerian stochastic fields: 32 fields converge to 27 % in the worst point;
  chemistry coordinate mapping (clustering) gives a 25× speed-up; within 14 %
  of experiment on Spray A; about 140 h for 1.5 ms (Pang et al. 2018 above).
- Codes: CMC, Lagrangian PDF (Ansys Fluent: Pei et al.), MMC (Sydney/Aachen).

### 2.5 Premixed-front models for spark ignition

- G-equation (level set of the flame front, Peters 2000), thickened flame
  (Colin et al., Phys. Fluids 2000), flame surface density (ECFM-LES).
- One comparison in an SI engine LES (Kazmouz et al., Sci. Tech. Energ.
  Transition 77 (2022) 10, [STET](https://www.stet-review.org/10.2516/stet/2022009))
  *(summary)*: the spread of peak pressure over cycles is 7.1 bar measured,
  8.16 bar with detailed chemistry per cell, 13.25 bar with thickened flame,
  3.63 bar with the G-equation (too small, kernel grows too fast).
- None is a clear winner; for us an SI engine is not in any plan yet.

### 2.6 Machine learning

Reviews of 2025–2026 (Shateri et al., [arXiv 2604.25617](https://arxiv.org/pdf/2604.25617);
Wu et al., [arXiv 2509.03347](https://arxiv.org/html/2509.03347v3)) list
neural-network tables (70 % less memory in one scramjet study), learned
chemistry integrators and learned sub-grid burning rates, all as research
demonstrations. Open problems: extrapolation outside the training data and
drift of mass and element conservation. Not production-ready; we watch it (it
fits the data-driven plan D205–D212 later, never as a default).

## 3. Ranking

| Family | Accuracy (evidence) | Cost (relative) | New fuel | Extinction, ignition | GPU fit |
|---|---|---|---|---|---|
| Tables FGM/FPV/UFPV | Spray A ID 5–20 %, LOL 1–14 %; Flame D good | 1 | rebuild table | weak on strong extinction; UFPV for ignition | very good (lookup) |
| Detailed chemistry per cell + clustering | Spray A ID and LOL too long (up to 67 %); fine in well-resolved LES | 10–40 (after 10× clustering) | none needed | yes, laminar-like | good (SF2b kernel) |
| EDC (standard) | furnace globals OK, MILD T +28 %, Flame D T +500 K, NO ×3 | 10–40 | none needed | poor in MILD | good |
| EDC local constants / PaSR | MILD T 9 %, furnace T 2.3–4.5 % | as EDC | none needed | better | good |
| Transported PDF (stochastic fields) / CMC | best LOL on Spray A (within 14 %) | 50–100+ | none needed | yes | fair (many fields) |
| G-equation / thickened flame / ECFM | engine-tuned | 1–5 | flame speed per fuel | separate models | good |
| ML tables or chemistry | research | <1 (if it works) | retrain | unknown | very good |

## 4. What we have and where it fits

| Have | State | Fits |
|---|---|---|
| SF1 species transport (#278) | merged | all models |
| SF2a 0D reactor, RODAS3, GRI 3.0, JL4, WD2 (#289) | merged | table generation (0D ignition, flamelets need a 1D solver added) |
| SF2b chemistry in flow cells, laminar transient, CPU = GPU bitwise (#298) | merged; refused with turbulence, sprays, steady, ranks (CHM003) | the "detailed chemistry per cell" family, already |
| SF3 EDC with finite-rate fine structures | being built (sibling thread) | furnace (D281: EDC **and** tables in the first version) |
| SF13 non-adiabatic flamelet / progress-variable tables | planned after SF2 | this survey's first recommendation |
| Spray A, Spray G, non-reacting; SP4/SP5 evaporation | merged | reacting Spray A next to tables |
| GPU, MPI ranks | species and variable density on ranks (#303) | tables are cell-local lookups: easiest GPU path of all |

Gaps the survey shows:
- No turbulence–chemistry closure for detailed chemistry per cell yet; SF3
  closes that with EDC; PaSR and local EDC constants are small additions on the
  same integrator.
- No clustering: every reacting cell is integrated (fine for tests, slow at
  millions of cells; Cuoci et al. 2024 and the LLNL multi-zone give 10–13×).
- No 1D flamelet solver (counterflow diffusion flame and freely propagating
  premixed flame with the SF2 kinetics); needed for any table. The SF2b flame
  speed case is already a 1D premixed flame in the 3D code.
- Preferential diffusion: SF2b has species Lewis numbers (needed for GRI 0.348
  m/s); tables for hydrogen need the H-radical enthalpy term (D281 asks for H₂
  to pure hydrogen).

## 5. Proposed decisions (cards, one at a time)

| No. | Question | Options | Recommendation |
|---|---|---|---|
| D360 | Which model next after SF3 EDC | (a) non-adiabatic flamelet tables (FGM premixed + FPV/UFPV diffusion), one table kernel for furnace and sprays; (b) chemistry clustering first; (c) transported PDF (stochastic fields); (d) a flame-surface model of the ECFM type | (a): best accuracy per core hour on Spray A, the furnace plan asks for tables in the first version (D281), cell-local lookups are the easiest GPU kernel, and it needs no per-fuel flame-speed tuning |
| D361 | Finite-rate upgrades on SF3 | (a) PaSR and Parente local EDC constants as keys, plus cell clustering of the SF2 integrator; (b) only clustering; (c) none for now | (a): the furnace data show standard EDC 7.4 % vs PaSR 2.3–4.5 %, MILD 28 % → 9 %; clustering 10–13× is needed by EDC, PaSR, well-mixed and table generation alike |
| D362 | Engine route | (a) reacting ECN Spray A (ID, lift-off, measured) with detailed chemistry per cell and with UFPV tables, side by side; (b) wait for SI engines; (c) none | (a): the open ECN data are the standard benchmark (D201), our non-reacting Spray A is merged, and the side-by-side reproduces the published ranking |
| D363 | High-fidelity reference model | (a) Eulerian stochastic fields after the tables; (b) CMC; (c) not planned | (a): closes chemistry exactly, re-uses the SF2 integrator per field, and on Spray A is the reference; later stage, not sized for v1.x |
| D364 | Hydrogen and ammonia in tables | (a) from the start: enthalpy cross-diffusion with the H radical, stretch as a table coordinate option, NO by its own transport equation; (b) later | (a): D281 asks for hydrogen blends to pure hydrogen; without it the table is 20 % fast on lean H₂ and adiabatic tables miss NO |
| D365 | Table storage and ML | (a) dense tables with tensor-train compression later for GPU memory, ML tables only as a research option (D205–D212); (b) neural-network tables now | (a): compression has an error bound; ML is not production-ready in any source |

Rough sizing (agent time, not sized in detail): 1D flamelet solver plus table
generator 1.5–2.5 h code (1–2 PRs); table lookup in flow cells with presumed
PDFs and the furnace and Flame D cases 2–3 h code plus 1–2 h runs; PaSR and
local EDC 30–60 min; clustering 1–2 h; reacting Spray A 2–3 h code plus 2–4 h
runs; stochastic fields 3–5 h code.

## 6. Validation targets (open data)

- ECN Spray A reacting: ignition delay and lift-off length, 750–1000 K
  ([ECN](https://ecn.sandia.gov)).
- TNF Sandia Flames D and F (piloted CH₄/air jet), already the SF3 target.
- Adelaide jet in hot coflow (HM1–HM3) for MILD and EDC constants.
- Lean H₂ premixed flames with published DNS (a-priori table checks).

## 7. Sources not read in full

The OSTI TFM paper (osti.gov purl 1579904) and the TRB record of a flamelet
study for hydrogen–ammonia SI engines did not load; two ScienceDirect papers
(transported PDF of n-dodecane, MMC on Spray A) were abstract-only. Nothing
from them is used above except their titles.
