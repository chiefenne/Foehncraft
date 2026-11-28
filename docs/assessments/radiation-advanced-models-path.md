# Advanced radiation models: upgrade path for accuracy and speed

Plan, 2026-11-28. Requested by the maintainer: "draft an update path
(accuracy and speed) for more advanced radiation models and card me on
that". Documentation only; **coding of each stage starts on the maintainer's
explicit go.** State read: `main` at 570b188, the
[radiation specification](../requirements/step29-radiation.md) (R1a done in
#254, #255, #257, #259, #281; R1b to R1g specified; D255 to D259, D269), D107,
D108, the [graphics algorithms assessment](graphics-algorithms-for-modules.md)
§4, the [steam furnace path](steam-furnace-development-path.md) (SF3 EDC
merged #351, SF4 = R1b + R1f + R1c, SF11 soot, SF16 solid fuels) and the
furnace rerun on the 2.1 M-cell layered mesh of 2026-11-25 (448 µs per bundle
on 16 threads against about 19 µs on the coarse 58 k-tetrahedron mesh on 4
threads: about 100 times more CPU per bundle, cause not profiled yet).

Clean room (D201): methods from papers and textbooks; references are exact
solutions, measurements, and line-by-line spectral results computed from a
spectroscopic database (the same accepted class as DNS, entered in D201 before
they appear in a table). Results of other CFD or radiation codes are not
reference numbers.

Stage ids use the prefix **RA** (radiation, advanced). Decision numbers
D366 to D372 are proposed (section 7, one card each); they are provisional
and checked against `main` again before this page is merged.

## 1. Answer in short

R1a gives an exact, reproducible, grey Monte Carlo kernel on CPU, ranks and
GPU. For an industrial furnace it lacks three things, in this order of
importance:

1. **Speed on real meshes.** The per-bundle cost grew about 100 times from
   the coarse to the layered furnace mesh. At 448 µs per bundle one call with
   one million bundles takes about 7.5 minutes of CPU on 16 threads; coupled
   to the flow every 20 to 50 outer iterations that is far more than the
   flow itself. This is the first thing to fix, before coupling (R1b), or the
   coupled furnace of SF5 is not affordable.
2. **Spectral gas properties.** Grey gas misses the band structure of H₂O
   and CO₂: errors of 20 to 50 % in the wall flux and the radiative source
   are typical (Modest & Mazumder 2021, chapter 20). WSGG with current
   coefficients brings that to about 5 to 20 %; full-spectrum k-distributions
   (FSCK) to a few percent of line-by-line. In Monte Carlo the spectral model
   costs almost nothing extra: each bundle samples one grey gas or one
   g-value, so the bundle count does not grow with the number of bands. This
   is the main structural advantage of Monte Carlo over discrete ordinates,
   whose cost grows linearly with the number of gases or quadrature points.
3. **Particles, soot, scattering and turbulence–radiation interaction** for
   sooting, oil and coal flames (D281) and for the radiant fraction of jet
   flames.

Recommended order (section 3): **RA0 speed foundation → RA1 coupling (R1b)
with cost control → RA2 WSGG → RA3 own line-by-line reference → RA4 FSCK /
SLW → RA5 scattering, soot and particles → RA6 non-grey and specular walls →
RA7 variance reduction → RA8 deterministic options (P1, FV-DOM) → RA9 GPU
wavefront and graphics benchmarks → RA10 turbulence–radiation interaction →
RA11 hybrid sphere steps (later).** RA0 to RA4 are the furnace's critical
path (SF4); RA5 is needed by SF11 and SF16.

## 2. What exists and what each furnace stage needs

| Capability | State | Needed by |
|---|---|---|
| Grey absorbing-emitting gas, grey diffuse walls, mirrors, open boundaries, track-length estimator, Sobol emission, integer tallies, ledger | done (R1a) | everything |
| Multi-rank hand-over, bitwise 1 to 4 ranks; GPU kernel with integer atomics | done (R1a, D258, D259, D269); real CPU = GPU run on a GPU still waits for the Colab check | SF9 package boiler |
| Paths of chosen bundles (`Domain::trace`) | done (#281) | gallery |
| Coupling to the energy equation, linearised emission, CHT wall flux | specified (R1b) | SF4, SF5, honeycomb |
| Scattering, next-event estimation | specified (R1c) | SF11, SF15, SF16 |
| Reciprocal estimator, MIS | specified (R1d) | thick or near-isothermal regions |
| Null-collision tracking on a majorant grid | specified (R1e) | speed benchmark |
| WSGG for H₂O/CO₂ | specified (R1f, Bordbar et al. 2014) | SF4 |
| FSCK, SLW, line-by-line, soot, particles, TRI, P1/DOM | not specified | SF5 accuracy, SF11, SF15, SF16, flame validation |

## 3. Ordered path

Gains are estimates from the cited literature and from the R1a numbers, not
measurements of Foehncraft; each stage measures its own. Code time is agent
time per pull request; runs separate.

| Stage | Content | Accuracy gain | Speed gain | Code / runs |
|---|---|---|---|---|
| **RA0** | **Speed foundation.** Profile the 2.1 M-cell furnace (where the 448 µs go: face tests of polyhedra, cache misses, hand-overs, roulette tails). Then: (1) per-cell face planes stored once (convex cells: exit = smallest positive plane distance, one division per face; the exact predicate of D225 only near ties), (2) cells and faces renumbered along a space-filling curve (Morton or Hilbert) so a bundle's walk touches neighbouring memory, (3) structure-of-arrays face data, (4) absorption tallies per bundle buffered and flushed per cell run. Run the pending Colab CPU = GPU check and record the GPU time per bundle. Same results to the quantum (the walk is unchanged in exact arithmetic; renumbering is a permutation of global ids that D93 already keys on) | none (bitwise or quantum-equal) | 5 to 20× on layered polyhedral meshes (estimate; the coarse mesh shows the attainable cost) | 1 to 2 PRs, 45 to 90 min; profiling runs 20 to 40 min on the owner's iMac |
| **RA1** | **R1b coupling with cost control.** As specified (Fleck–Cummings linearisation, running average, CHT wall flux, ledger rows), plus: bundle count chosen from a target statistical error of the wall flux and the source (D69 health number) instead of a fixed number; call interval adapted to the change of the temperature field between calls (no call while the source is converged within its own noise); the source frozen between calls with the linearised emission keeping the energy equation stable | coupled temperatures and wall fluxes (the furnace needs it) | radiation share of a steady furnace run held at 10 to 20 % of the flow time (target) | 2 PRs as specified, plus 1 for the controller |
| **RA2** | **WSGG (R1f) with current coefficient sets.** One grey gas sampled per bundle with probability a_j(T) of the emitter; absorption per cell with κ_j(T, x_H₂O, x_CO₂). Coefficient sets: Bordbar et al. 2014 (H₂O/CO₂ any ratio, line-by-line fitted), Bordbar et al. 2020/2021 (with soot; citation to check), Ehlmé et al. 2025 (to 5000 K, Chalmers). Keyed by `[radiation.gas] model = "wsgg"` with `coefficients = "..."` and provenance in the run record | grey → WSGG: wall flux error typically from 20 to 50 % down to 5 to 20 % against line-by-line in non-homogeneous media; the 2024 comparison of new WSGG sets (arXiv 2411.18467) shows the spread between sets in oxy-fuel | none per bundle; same bundle count | 1 PR, 45 to 60 min |
| **RA3** | **Own line-by-line reference (offline tool).** Absorption coefficients from HITEMP 2010 (H₂O, CO₂) and HITRAN (other species) with Voigt profiles on a fine wavenumber grid; spectral Monte Carlo (wavenumber sampled per bundle from the emission spectrum) in 1D slabs and 3D boxes. Two uses: the reference for RA2/RA4 and the generator of our own WSGG and FSCK tables, so the solver depends on no external fit | the reference itself (the accepted class of D201) | n/a (offline, hours per table) | 2 PRs, 60 to 90 min each; table generation runs 1 to 2 h |
| **RA4** | **FSCK and SLW.** Full-spectrum correlated-k (Modest & Zhang 2002; Modest 2003) with the multi-scale form for strongly non-homogeneous flames, tables from RA3 (or the look-up approach of Wang, Ge, Modest & He 2016 rebuilt from our own data); SLW (Denison & Webb 1993) and its one-gas form SLW-1 (Solovjov et al. 2017) as the cheap alternative. In Monte Carlo: one g-value per bundle (Wang & Modest 2007) | within about 2 to 5 % of line-by-line in non-homogeneous furnaces | table look-ups only; the bundle count is set by the noise, not by the spectrum | 2 PRs, 60 to 90 min each |
| **RA5** | **Scattering, soot, particles (R1c plus).** Isotropic and Henyey–Greenstein scattering, next-event estimation towards walls with ratio tracking (Novák et al. 2014); soot absorption in the Rayleigh limit with κ ∝ f_v T (refractive index Chang & Charalampopoulos 1990) inside WSGG/FSCK (Bordbar 2020 or own RA3 tables with soot); parcels of sprays and coal (SP layer) as absorbing and scattering particles with efficiencies from our own Mie routine (Bohren & Huffman 1983 algorithm) per size class; ash later (SF17) | sooting and particle-laden flames (oil SF15, coal SF16) become possible at all; soot can double the radiant flux of an oil flame | next-event estimation lowers the wall-flux noise by 2 to 10× in scattering media (estimate) | 2 to 3 PRs |
| **RA6** | **Non-grey and specular walls.** Emissivity per WSGG gas or per FSCK g-band (two or three bands), specular plus diffuse reflection by a fraction, semi-transparent windows (inspection ports, quartz walls of lab burners) | refractory and fouled tube surfaces (emissivity changes from 0.6 to 0.9 with deposits, SF17) | none | 1 PR, 45 min |
| **RA7** | **Variance reduction (R1d plus).** Reciprocal estimator (Tessé et al. 2002) for net exchange in thick and near-isothermal regions; MIS between forward and reciprocal estimators (Veach & Guibas 1995); importance in the emission (more bundles where the wall flux is wanted, with weights); spatial filtering of the source only for the coupling (conservative, the ledger stays exact). Graphics benchmarks: path guiding (Müller et al. 2017) for the directional sampling towards tubes | noise at equal cost: 3 to 30× lower in thick regions and in wall flux maps (estimate from Tessé 2002, Zhang & Modest 2003 for reciprocal MC) | the same error with 3 to 30× fewer bundles | 2 PRs |
| **RA8** | **Deterministic options.** P1 on the existing scalar matrix path (robust in thick, sooty or coal flames; wrong in thin gas, so not the default); FV-DOM (Raithby & Chui 1990; Chui & Raithby 1993) with S_N and FT_N quadratures on unstructured meshes, ray-effect mitigation by the random ordinate method (Li, Tang & Yang 2024) or random source iteration (Fu, Li & Tang 2025), WSGG/FSCK per gas. MC stays the default; DOM and P1 are speed options and cross-checks against the MC result, and the users of industrial codes expect them | P1: poor in thin media; FV-DOM: ray effects and false scattering, mitigated by randomised ordinates; within a few % of MC with fine angles | P1 is 10 to 100× cheaper than MC per call; DOM per call costs (angles × gases) scalar sweeps, e.g. 80 × 5 = 400 for S8 with 4+1 WSGG gases | 2 to 3 PRs |
| **RA9** | **GPU wavefront and graphics benchmarks.** Wavefront scheduling (Laine, Karras & Aila 2013): bundles sorted by cell between steps so threads of a warp walk the same cells; persistent threads; null-collision tracking on a majorant grid (R1e, Woodcock et al. 1965; Galtier et al. 2013) benchmarked against the walk; BVH ray casting only for transparent regions (the heat-exchanger case: surface-to-surface through air skips the cell walk). RT cores stay out (float64, D107 (6)) | none (same estimator) | GPU kernels of this type report 10 to 50× over a multi-core CPU (Silvestri & Pecnik 2018 for DNS-coupled MC; PUMI-Tally 2025 about 20× on A100/H100 for unstructured track-length tallies) | 2 PRs; GPU runs on Colab or the cloud |
| **RA10** | **Turbulence–radiation interaction.** RANS: emission TRI from a transported temperature variance and the self-correlation of the Planck function (Li & Modest 2003; Coelho 2007 review), absorption TRI neglected first (optically thin fluctuation approximation, valid for most furnaces); with the stochastic fields of D361 or a PDF method later: emission and absorption per field or particle, which is exact for emission TRI in Monte Carlo (Wang & Modest 2007). LES (L5): resolved fluctuations in the instantaneous fields, subfilter part neglected (Roy et al.; Coelho 2012) | radiant fraction of jet flames: TRI raises emission by 30 to 50 % (Flame D); in large furnaces 10 to 20 % on the wall flux | none (RANS correction is a field) | 1 to 2 PRs |
| **RA11** | **Hybrid sphere steps** in optically thick regions (L15 (2), Fleck–Canfield 1984 lineage and the 2022 to 2024 walk-on-spheres estimators): published-research work, after RA7 shows where the time goes | equal to MC | potentially large in thick slag, ash or porous solids | not sized yet |

Dependencies: RA0 → RA1 → (RA2 → RA4), RA3 feeds RA2 and RA4 references;
RA5 needs RA1 and SP1; RA7 after RA1; RA8 independent after RA2; RA9 after
RA0 and the Colab check; RA10 after RA2 and SF3; RA11 after RA7.

## 4. Coupling cost per outer iteration

Per radiation call, Monte Carlo costs N × t_b / threads; R1a measured
t_b ≈ 19 µs (coarse furnace) and 448 µs (layered 2.1 M-cell furnace). A
steady furnace run of a few thousand outer iterations at about 1 s each needs
the radiation share near 10 to 20 % of that. With RA0 bringing t_b to about
25 to 50 µs on the fine mesh, 1 M bundles on 16 threads cost 2 to 3 s per
call; called every 20 outer iterations that is 10 to 15 % of the flow time,
and the adaptive interval of RA1 calls less as the field settles. On the GPU
(RA9) the call drops below the flow iteration.

FV-DOM costs per call one sweep per angle and per gas: S8 has 80 angles, a
4-gas WSGG adds a transparent window: 400 sweeps of a scalar equation, each
comparable to one scalar solve. P1 costs one scalar solve per gas. This is
why DOM is the speed option only with grey or few-gas models, and why FSCK
in DOM needs 8 to 16 quadrature points (×80 angles), while Monte Carlo needs
none extra.

## 5. Validation cases with open data

| Stage | Case | Data | Status for D201 |
|---|---|---|---|
| RA0 | R1a cases 1 to 6, the furnace showcase | exact / own quadrature; result to the quantum | in |
| RA1 | Honeycomb quarter channel (R1b); slab with conduction and radiation (Viskanta & Grosh 1962 type, exact in the thin and thick limits) | exact limits, own resolved reference | in |
| RA2, RA4 | H₂O/CO₂ slabs and the non-isothermal, non-homogeneous slab cases of Bordbar et al. 2014 (line-by-line tables in the paper); 3D enclosures of Goutière, Liu & Charette 2000 and Coelho 2002 with statistical narrow-band references | published line-by-line or narrow-band results, and our own RA3 line-by-line | needs D201 entries (spectral references, same class as DNS) |
| RA3 | HITEMP 2010 (Rothman et al. 2010), HITRAN | spectroscopic database, free with registration | licence check before download; data not redistributed in the repo |
| RA5 | 3D cube with isotropic and anisotropic scattering (Crosbie & Schrenker 1984; Hsu & Farmer 1997); sooting ethylene jet flames of the International Sooting Flame workshop (soot volume fraction and radiant flux) | exact or published measurements | ISF data needs its terms checked |
| RA6 | Enclosures with specular walls (Howell catalogue, exact view factors with mirrors) | exact | in |
| RA8 | Same cases as MC; ray-effect case: small hot square in a cold enclosure (Chai, Lee & Patankar 1993) | exact / MC reference | in |
| RA10 | Sandia Flame D radiant fraction and radiative flux profiles (TNF workshop, Barlow; Frank et al. 2000) | measurements, open on the TNF site | needs the TNF terms checked; the owner already asked for Flame D data for SF3 |
| SF5 furnace | Chalmers 100 kW oxy-fuel and air propane unit: in-flame temperature, gas composition and radiative intensity profiles (Andersson et al. 2008, 2011); IFRF natural gas furnaces (Weber et al.) with wall flux | open in journal papers | D201 entries; the pilot furnace of SF5 |

## 6. Graphics-method benchmarks (owner rule 2026-11-21)

Every stage names the computer-graphics method it is measured against:
RA0 face-plane walk against the predicate walk and against a BVH ray cast
through the cell hierarchy of D108; RA7 path guiding and MIS against the
forward estimator; RA9 wavefront scheduling against the megakernel, and
null-collision tracking against the walk (R1e); denoisers of real-time
rendering (spatiotemporal filters) are measured as post-processing of wall
flux maps only, never inside the coupled source, because they do not
conserve energy. RT-core traversal is excluded by D107 (6) (float64); a
float32 traversal stays the later item P9.

## 7. Questions (one card each, in this order)

| # | Question | Options | Recommended | Answer |
|---|---|---|---|---|
| Q1 (D366) | What comes first: speed or accuracy? | (a) RA0 speed foundation first, then R1b; (b) R1b coupling first as specified, speed after; (c) spectral models (RA2) first | (a): at 448 µs per bundle the coupled furnace is not affordable, and RA0 changes no result | **Answered 2026-11-28 (D366): (a) speed first.** RA0 (profiling, face-plane walk, renumbering, Colab GPU check) before R1b |
| Q2 (D367) | Spectral model ladder | (a) WSGG default now (RA2), FSCK/SLW as the accurate option (RA4); (b) FSCK as default from the start; (c) WSGG only | (a): WSGG is enough for SF4 and cheap to verify; FSCK follows once RA3 gives our own tables | open |
| Q3 (D368) | Spectral reference data | (a) own line-by-line tool from HITEMP, own WSGG/FSCK tables; (b) published line-by-line results and coefficient sets only; (c) both | (c): published sets give the first WSGG at once; our own tables remove the dependence on external fits | open |
| Q4 (D369) | Deterministic solvers | (a) P1 and FV-DOM with random-ordinate ray-effect fix as options, MC default; (b) MC only; (c) DOM default for speed | (a): users expect them, they cross-check MC, and P1 is the cheap option for thick sooty and coal flames | open |
| Q5 (D370) | Turbulence–radiation interaction | (a) RANS emission TRI from temperature variance first, full TRI with stochastic fields or PDF later; (b) full PDF Monte Carlo TRI from the start; (c) none | (a): most of the effect for little code; the full form comes with D361 stochastic fields | open |
| Q6 (D371) | Lead validation set | (a) Bordbar/Goutière/Coelho non-grey benchmarks, Sandia Flame D radiant fraction, Chalmers 100 kW unit; (b) add IFRF furnaces now; (c) analytic cases only until SF5 | (a): open data that covers spectra, TRI and a real furnace | open |
| Q7 (D372) | Scope of the GPU and graphics benchmarks | (a) wavefront, cell sorting, null-collision, path guiding; RT cores stay out (float64); (b) also a float32 RT-core traversal option (amends D107 (6), brings P9 forward) | (a): float64 and bitwise tallies are design rules; RT cores would break both | open |

## 8. Sizing

RA0 to RA4 (the furnace's radiation): about 8 to 11 pull requests of 45 to
90 minutes of agent code each, 6 to 12 h of code in all; profiling, table
generation and validation runs 4 to 8 h, mostly as background runs. RA5 to
RA10: about 10 to 13 pull requests, 9 to 15 h of code; validation runs of
flames 4 to 10 h. RA11: not sized yet.

## 9. References

- Bohren, C. F. & Huffman, D. R. (1983) *Absorption and Scattering of Light by Small Particles*. Wiley.
- Bordbar, M. H., Węcel, G. & Hyppänen, T. (2014) A line by line based weighted sum of gray gases model for inhomogeneous CO₂–H₂O mixture in oxy-fired combustion. *Combust. Flame* 161.
- Bordbar, H., Hostikka, S. and co-authors (2020 to 2021) WSGG models for CO₂–H₂O mixtures with soot and at other pressures. *Int. J. Heat Mass Transfer* / *Combust. Flame*. Exact citation to be checked against the paper before RA2 is built.
- Chai, J. C., Lee, H. S. & Patankar, S. V. (1993) Ray effect and false scattering in the discrete ordinates method. *Numer. Heat Transfer B* 24.
- Chang, H. & Charalampopoulos, T. T. (1990) Determination of the wavelength dependence of refractive indices of flame soot. *Proc. R. Soc. Lond. A* 430.
- Chui, E. H. & Raithby, G. D. (1993) Computation of radiant heat transfer on a nonorthogonal mesh using the finite-volume method. *Numer. Heat Transfer B* 23.
- Coelho, P. J. (2002) Numerical simulation of radiative heat transfer from non-gray gases in three-dimensional enclosures. *JQSRT* 74.
- Coelho, P. J. (2007) Numerical simulation of the interaction between turbulence and radiation in reactive flows. *Prog. Energy Combust. Sci.* 33.
- Coelho, P. J. (2012) Turbulence–radiation interaction: from theory to application in numerical simulations. *J. Heat Transfer* 134.
- Crosbie, A. L. & Schrenker, R. G. (1984) Radiative transfer in a two-dimensional rectangular medium exposed to diffuse radiation. *JQSRT* 31.
- Denison, M. K. & Webb, B. W. (1993) A spectral line-based weighted-sum-of-gray-gases model for arbitrary RTE solvers. *J. Heat Transfer* 115.
- Ehlmé, E., Gunnarsson, A., Normann, F. & Andersson, K. (2025) Updated weighted-sum-of-gray-gases model parameters for a wide range of water and carbon dioxide concentrations and temperatures up to 5000 K. *ACS Omega* 10(3), 2978–2985.
- Frank, J. H., Barlow, R. S. & Lundquist, C. (2000) Radiation and nitric oxide formation in turbulent non-premixed jet flames. *Proc. Combust. Inst.* 28.
- Fu, J., Li, L. & Tang, M. (2025) Random source iteration method: mitigating the ray effect in the discrete ordinates method. arXiv:2502.15454.
- Galtier, M. et al. (2013) Integral formulation of null-collision Monte Carlo algorithms. *JQSRT* 125.
- Goutière, V., Liu, F. & Charette, A. (2000) An assessment of real-gas modelling in 2D enclosures. *JQSRT* 64.
- Hasan, F. et al. (2025) GPU acceleration of Monte Carlo tallies on unstructured meshes in OpenMC with PUMI-Tally. arXiv:2504.19048.
- Hsu, P.-F. & Farmer, J. T. (1997) Benchmark solutions of radiative heat transfer within nonhomogeneous participating media using the Monte Carlo and YIX method. *J. Heat Transfer* 119.
- Laine, S., Karras, T. & Aila, T. (2013) Megakernels considered harmful: wavefront path tracing on GPUs. *High-Performance Graphics*.
- Li, G. & Modest, M. F. (2003) Importance of turbulence–radiation interactions in turbulent diffusion jet flames. *J. Heat Transfer* 125.
- Li, L., Tang, M. & Yang, Y. (2024) Random ordinate method for mitigating the ray effect in radiative transport equation simulations. arXiv:2407.12527.
- Modest, M. F. (2003) Narrow-band and full-spectrum k-distributions for radiative heat transfer: correlated-k vs. scaling approximation. *JQSRT* 76.
- Modest, M. F. & Mazumder, S. (2021) *Radiative Heat Transfer*, 4th ed. Academic Press.
- Modest, M. F. & Zhang, H. (2002) The full-spectrum correlated-k distribution for thermal radiation from molecular gas–particulate mixtures. *J. Heat Transfer* 124.
- Müller, T., Gross, M. & Novák, J. (2017) Practical path guiding for efficient light-transport simulation. *Comput. Graph. Forum* 36(4).
- Novák, J., Selle, A. & Jarosz, W. (2014) Residual ratio tracking for estimating attenuation in participating media. *ACM TOG* 33(6).
- Raithby, G. D. & Chui, E. H. (1990) A finite-volume method for predicting radiant heat transfer in enclosures with participating media. *J. Heat Transfer* 112.
- Rothman, L. S. et al. (2010) HITEMP, the high-temperature molecular spectroscopic database. *JQSRT* 111.
- Silvestri, S. & Pecnik, R. (2018) A fast GPU Monte Carlo radiative heat transfer implementation for coupling with direct numerical simulation. arXiv:1810.00188 (*J. Comput. Phys. X* 2019).
- Solovjov, V. P., Andre, F., Lemonnier, D. & Webb, B. W. (2017) The rank correlated SLW model of gas radiation in non-uniform media. *JQSRT* 197.
- Tessé, L., Dupoirieux, F., Zamuner, B. & Taine, J. (2002) Radiative transfer in real gases using reciprocal and forward Monte Carlo methods and a correlated-k approach. *Int. J. Heat Mass Transfer* 45.
- Veach, E. & Guibas, L. J. (1995) Optimally combining sampling techniques for Monte Carlo rendering. *SIGGRAPH*.
- Wang, A. & Modest, M. F. (2007) Spectral Monte Carlo models for nongray radiation analyses in inhomogeneous participating media. *Int. J. Heat Mass Transfer* 50.
- Wang, C., Ge, W., Modest, M. F. & He, B. (2016) A full-spectrum k-distribution look-up table for radiative transfer in nonhomogeneous gaseous media. *JQSRT* 168.
- Woodcock, E. et al. (1965) Techniques used in the GEM code for Monte Carlo neutronics calculations in reactors and other systems of complex geometry. ANL-7050.
- Comparison of new WSGG sets for carbon-capture conditions (2024) arXiv:2411.18467.
