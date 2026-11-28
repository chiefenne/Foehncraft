# General multiphase: any number of phases and several bubble sizes per phase (D373 to D376)

Status: assessment 2026-11-28; D373 to D376 answered the same day (section 8); **no coding before the maintainer's go**. Assessment for the later multiphase work (L3). Clean-room rules of D201 apply: only open papers are cited, the code of GPL projects (OpenFOAM, the HZDR addon, OpenQBMM) was not read, and every equation below names its source with a section or equation number. What could not be read is listed in section 10.

## 0. The request

Maintainer, 2026-11-28 (tube boiling thread): "for later multi-phase add to the reqs we want arbitrary number of phases (although practical use will limit itself to few) so a general implementation. if several bubble sizes are possible (also per phase) that would be good. a researcher i follow is Alessio Roccon, check his and related research what is leading edge for the request i put here." (spelling corrected)

So two requirements: (a) N phases with no fixed limit in the code; (b) several bubble or drop sizes, also within one phase. Plus: where does the leading edge stand, and how does this sit on the MP3 two-fluid core.

## 1. Short answer

1. **Roccon's work is the leading edge of interface-resolved DNS of drops and bubbles, not of N-phase models.** His codes (FLOW36, pseudo-spectral Cahn-Hilliard; MHIT36, finite-difference conservative Allen-Cahn on GPUs, MIT licence) resolve every interface of two fluids. The number of drops is arbitrary and the size distribution comes out of the simulation by counting drops. What we take from it is physics: the drop size distribution splits at the Kolmogorov-Hinze scale into a coalescence range (d^-3/2) and a breakup range (d^-10/3), and DNS of this kind is now the source of breakup and coalescence rates for the models an industrial code needs (section 3).
2. **For an industrial solver the general form is a multi-fluid Euler-Euler model with N phases, where a phase is a velocity field, and a size model per dispersed phase.** This is how the inhomogeneous MUSIG model is built (N velocity groups, M size groups each; Frank et al. 2005, Eqs. 8 to 12; Krepper et al. 2007, Eqs. 5 to 9), and how the moment methods with a size-conditioned velocity are built (Heylmun, Fox and Passalacqua 2021; Li et al. 2019). Several bubble sizes "also per phase" are then natural: size groups inside a phase share its velocity, and sizes that must move differently (the lift force changes sign near 5.8 mm for air in water, Krepper et al. 2007 Sec. 3.1) go to another phase (section 5).
3. **MP3 fits without change.** MP3 builds one momentum equation per phase, drag by partial elimination and one pressure equation from the volume continuity of both phases (D352 to D356). Each of these generalises to N phases directly: the partial elimination becomes a small N×N solve per cell, the pressure equation sums over N phases (section 4.2). I recommend keeping MP3 and MP4 as decided and making the N-phase core the next stage (MP5), with the two-phase rows bit for bit as the case N = 2, then size groups (MP6), moments later (MP7), and the bridge between resolved and dispersed interfaces after that (MP8). D373 to D376 confirm this order.

## 2. What exists today

| Piece | State | Relevant to N phases |
|---|---|---|
| MP1 VOF (D293, #317) | merged | two phases, one α; case `[[vof.phase]]` is already an array but `validate` takes exactly two |
| MP2 phase change in VOF (D351, #358) | merged 2026-11-28 | one liquid, one vapour |
| MP3 two-fluid core (D352 to D356) | in progress, PR not up yet | two momentum equations, partial elimination drag, one pressure; one bubble size; files `foehncraft-flow/src/twofluid.rs`, `foehncraft-kernels/src/twofluid.rs` (to be re-read when the PR opens) |
| MP4 wall boiling (Kurul and Podowski) | next | bubble departure size and nucleation sites feed the size model later |
| L12 parcels (sprays) | merged, many stages | a Lagrangian route to any size distribution |
| D194 several fluids and materials | MF-A merged | per-cell region index and properties as fields; a template for "no fixed limit in the code" |
| SF7b circuits (1D) | merged | drift-flux void, not a size model |

D293 (3), the conservative momentum time term with mass fluxes ρ_f F_f from the advected α, is the two-phase case of the "consistency of mass and momentum transport" that Huang, Lin and Ardekani (2022, Sec. 2.2 and Eq. 12) require for N phases. Our VOF is therefore already on the right principle.

## 3. Alessio Roccon and related groups

### 3.1 Methods

| Code / paper | Method | Phases | GPU | Licence |
|---|---|---|---|---|
| Roccon, "A GPU-ready pseudo-spectral method for direct numerical simulations of multiphase turbulence", Procedia Computer Science (EuroHPC user day), 2024, arXiv:2406.01322 (read) | FLOW36: Navier-Stokes plus Cahn-Hilliard, pseudo-spectral (Fourier, Chebyshev). Sec. 2.1.1: ∂φ/∂t + u·∇φ = (1/Pe) ∇²μ (Eq. 1), μ = φ³ − φ − Ch²∇²φ (Eq. 4), surface tension as the divergence of the Korteweg stress (3/√8)(Ch/We)∇·τ_c, τ_c = \|∇φ\|²I − ∇φ⊗∇φ (Eq. 8) | two (matched density and viscosity in this paper, Sec. 2.1.2) | MPI plus OpenACC and cuFFT, about 16× on Leonardo against all CPU cores | article CC BY-NC-ND; code in-house |
| Mangani, Soligo, Roccon, Soldati, "Influence of density and viscosity on deformation, breakage, and coalescence of bubbles in turbulence", arXiv:2309.13995 (read, venue not confirmed) | same code with density and viscosity contrast: ρ(φ) and η(φ) linear in φ (Eqs. 10, 11); profile-corrected flux f_p = (λ/Pe)[∇²φ − (1/(√2 Ch)) ∇·((1 − φ²)∇φ/\|∇φ\|)] (Eq. 7) | two (bubbles: ρ_r down to 0.001, η_r 0.01 to 100) | as above | arXiv |
| Roccon, Enzenberger, Zaza, Soldati, "MHIT36: A phasefield code for GPU simulations of multiphase homogeneous isotropic turbulence", Comput. Phys. Commun. 316 (2025) 109804 (first half read) | second-order finite differences, staggered, projection with FFT Poisson solve; conservative Allen-Cahn of the accurate conservative diffuse-interface type: ∂φ/∂t + ∇·(uφ) = ∇·[Γ(ε∇φ − ¼(1 − tanh²(ψ/2ε)) ∇ψ/\|∇ψ\|)] (Eq. 9), ψ = ε ln((φ + ξ)/(1 − φ + ξ)) (Eq. 10), f_σ = μ_φ∇φ (Eq. 3) | two | OpenACC, cuFFT, cuDecomp; scaling to 1024 GPUs stated in the abstract | MIT, github.com/MultiphaseFlowLab/MHIT36 |

Reviews by the group (metadata read, full text not reachable): Soligo, Roccon, Soldati, "Turbulent flows with drops and bubbles: what numerical simulations can tell us" (Freeman Scholar Lecture), J. Fluids Eng. 143 (2021) 080801; Roccon, Zonta, Soldati, "Phase-field modeling of complex interface dynamics in drop-laden turbulence", Phys. Rev. Fluids 8 (2023) 090501 (review, breakup, coalescence, size distribution, surfactants, heat and mass transfer, per its abstract).

### 3.2 Physics that matters for a size model

- **Size distribution and the Hinze scale.** Mangani, Roccon, Zonta and Soldati, "Heat transfer in drop-laden turbulence", J. Fluid Mech. (2024), read as arXiv:2309.13976v1, Sec. 3.2: the equivalent diameter d_eq = (6V/π)^(1/3) (Eq. 3.1) is counted per drop; the Kolmogorov-Hinze scale d_H+ = 0.725 (We/Re_τ)^(-3/5) \|ε_c\|^(-2/5) (Eq. 3.2) separates a coalescence-dominated range with d^-3/2 from a breakup-dominated range with d^-10/3 (Fig. 4). The same two laws are cited there to Soligo et al. (J. Fluid Mech. 881, 2019), which I could not read; that is the second source in name only.
- **Bubbles against drops.** Mangani et al. (arXiv:2309.13995, Sec. IV, verbatim part): the density ratio changes the number of bubbles little, the viscosity ratio a lot (higher η_r, fewer breakups). So a breakup kernel for bubbles cannot be taken from drops without the viscosity ratio.
- **What this means for us.** These are two-fluid DNS at Re_τ = 300 on 512 × 256 × 513 points (Mangani et al. 2023, Table I). They give us test laws (the slopes around d_H) for a homogeneous check of any size model, and they are where breakup and coalescence rates come from today. They are not a model we can run on a boiler tube.

### 3.3 Related groups at the leading edge (2024 to 2026)

- **Breakup rates from DNS.** Morón, Cannon, Vela-Martín, Avila, "A causal model of drop breakup in turbulence", arXiv:2605.07504v3 (2026, read, abstract and Sec. 4): breakup at low Weber number is caused by outer eddies that stretch the interface faster than capillarity can return the energy; a breakup-rate equation (Eq. 4.9) fits drops and bubbles over a wide range of We with one constant (κ̂_∞ ≈ 10). A candidate kernel for size groups, and a modern alternative to Luo and Svendsen (1996).
- **Small fragments.** Crialesi-Esposito, Boffetta, Brandt, Chibbaro, Musacchio, "How small droplets form in turbulent multiphase flows", arXiv:2402.19121 (2024, abstract read): small drops come from the internal dynamics of a breaking large drop, not from strong eddies. This is the daughter-size distribution a kernel needs.
- **Many small bubbles.** de Wit, Adelerhof, Freitas, Kunnen, Clercx, Toschi, "Dynamics of small bubbles in turbulence in non-dilute conditions", arXiv:2510.03968v2 (2026, abstract read): with four-way coupling, excluded volume between bubbles changes their statistics from about 1 % volume fraction.
- **Moment methods on GPUs, with heat and mass transfer.** Radhakrishnan and Bryngelson, "A sub-grid-scale model for polydisperse bubbly flows with heat and mass transfer", arXiv:2607.27426 (2026, abstract read): conditional hyperbolic moments in MFC, bubble pressure and vapour mass per node, 1.5 % rms against 40 Euler-Lagrange realisations. Vaca-Revelo et al., arXiv:2511.21031 (2025): ensemble-averaged against volume-averaged bubbly models on GPUs, 16× on four A100 against 64 cores. Both are compressible cavitation codes; the method family transfers, the solver does not (same verdict as our MFC assessment).
- **Moments, state of the art.** Denner, Campos, Blais, "QMOMKit: A C++ library for quadrature-based moment methods", arXiv:2610.36189 (2026, Sec. 1 read): surveys QMOM, DQMOM, EQMOM, HyQMOM, GQMOM and CQMOM with references (section 5.4 below); Apache 2.0. Posey, Fox, Houim, arXiv:2603.13697 (2026): GQMOM for continuous size distributions in high-speed flows.
- **Industrial polydisperse two-fluid (HZDR).** The inhomogeneous MUSIG papers (section 5.3); Lehnigk et al., "An open-source population balance modeling framework for the simulation of polydisperse multiphase flows", AIChE J. (2022) e17539 (not readable from here); the morphology-adaptive multifield model: Schlegel, Meller, Krull, Lehnigk, Tekavčič, "OpenFOAM-Hybrid: a morphology adaptive multifield two-fluid model", Nucl. Sci. Eng. (2022), doi 10.1080/00295639.2022.2120316 (project page read, paper not read; code GPL v3, not read).

## 4. An arbitrary number of phases

### 4.1 Interface-resolving with N phases

| Approach | How N phases are carried | Weak points | Sources |
|---|---|---|---|
| Multi-material VOF / moment of fluid | one α per material, interface reconstruction per cell | reconstruction "becomes more and more complicated as the number of phases increases" | Schofield et al. JCP 228 (2009) 731; Schofield et al. IJNMF 63 (2010) 931; Pathak and Raessi JCP 307 (2016) 550; all as cited and summarised by Huang et al. 2022, Sec. 1 |
| N level sets | one level set per phase | overlaps and voids after independent advection, mass loss | Huang et al. 2022, Sec. 1 |
| N-phase phase field | N order parameters, Σ C_p = 1 | needs reduction consistency (no fictitious phases appear) and consistent mass and momentum transport | Boyer and Minjeaud, Math. Models Methods Appl. Sci. 24 (2014) 2885; Dong, JCP 361 (2018) 1; Huang, Lin, Ardekani, J. Comput. Appl. Math. 408 (2022) 114116 (arXiv:2010.01099, read) |

Huang et al. (2022) is the cleanest N-phase statement: phase-field equations ∂φ_p/∂t + ∇·(uφ_p) = ∇·J_p with Σ_p J_p = 0 (Eq. 7), so Σ φ_p stays at 2 − N and ∇·u = 0 follows (Eq. 8); the consistent mass flux m = Σ_p (ρ_p/2)(u + uφ_p − J_p) (Eq. 12) carries the momentum, with the inertia in conservative form (Eq. 15); pairwise surface tensions σ_pq. The companion paper (Huang et al., arXiv:2101.04252, read in part) adds M dissolved components and gives the N-phase mobility and chemical potential with pairwise mixing coefficients λ_pq = (3/(2√2)) η σ_pq (Eqs. 6 to 10). The scheme keeps Σ = 1 in every cell for any N (abstract). The largest case I could verify is three phases (three-phase Rayleigh-Taylor, Sec. 5.2.1).

**Computer-graphics candidates (maintainer's benchmark rule).** Graphics has solved "many liquids" for animation: Losasso, Shinar, Selle and Fedkiw, "Multiple interacting liquids", ACM TOG 25 (2006) (read): one particle level set per region, advected independently, then overlaps and voids removed "by subtracting the average of the smallest two φ_i from all of the φ_i" (Sec. 3.1), density and viscosity jumps by modified coefficients in the pressure Poisson equation (Sec. 4.1), up to 10 phases shown (Fig. 7). This projection is cheap and local, a direct candidate to benchmark against N α fields with renormalisation. Da, Batty and Grinspun, "Multimaterial mesh-based surface tracking", ACM TOG 33 (2014), a non-manifold triangle mesh with a material label on each side (title and listing verified, page not readable); Ren, Li, Yan, Lin, Bonet and Hu, "Multiple-fluid SPH simulation using a mixture model", ACM TOG (2014), up to eight fluids with volume fractions and drift velocities (Table II), which is the Euler-Euler mixture model in particle form. Keenan Crane's 2020 to 2026 papers (publication list read) have nothing on multi-fluid tracking; "Subgrid marching tetrahedra" (TOG 2026, several surface sheets per cell) and "Winding numbers on discrete surfaces" (TOG 2023, region labels) are candidates for extracting and labelling several interfaces in output and for counting bubbles. A benchmark of these against multi-material VOF and the N-phase phase field belongs to the interface-resolved stage, not to MP5.

### 4.2 Euler-Euler with N phases (the industrial route)

Equations per phase k = 1..N (the MP3 equations with an index): continuity ∂(α_k ρ_k)/∂t + ∇·(α_k ρ_k u_k) = Σ_l Γ_lk; momentum with −α_k∇p and the interfacial forces of every pair (k, l); energy per phase when needed; Σ_k α_k = 1. Frank et al. (2005, Eqs. 1 to 3) and Krepper et al. (2007, Eqs. 5 to 7) write exactly this for one liquid and N gas velocity groups. Four design points:

1. **Implicit drag for N phases (my derivation, to be checked against the MP3 code).** After assembly, each cell and velocity component holds (a_k + Σ_{l≠k} K_kl) u_k − Σ_{l≠k} K_kl u_l = H_k − α_k V ∇p, with K_kl = K_lk ≥ 0 the drag coefficients of the pairs. The matrix M = diag(a) + (weighted graph Laplacian of K) is symmetric with a positive diagonal and diagonally dominant, so it is positive definite and a Cholesky solve of size N per cell is safe. This is the two-phase partial elimination of MP3 extended to N; Hanimann et al. (2021, below) describe the two-phase PEA and the Np × Np per-cell solve of SINCE (Karema and Lo, 1999), from a summary I could not verify. When drag only couples each dispersed phase to the continuous one (the usual bubbly case), M is an arrowhead and the solve costs O(N). The same M⁻¹ must enter the face (Rhie-Chow) velocities, otherwise the result depends on the time step: that is the subject of Hanimann, Mangani, Darwish, Casartelli, Vogt, "A consistent and implicit Rhie-Chow interpolation for drag forces in coupled multiphase solvers", Int. J. Turbomach. Propuls. Power 6 (2021) 7, CC BY (title and abstract read). It matters for our coupled solver (step 22) as well.
2. **One pressure from volume continuity.** Σ_k ∇·(α_k u_k) = Σ_k (Γ terms)/ρ_k, as in MP3 and in D293 (2); with N phases the face coefficient becomes Σ_k Σ_l α_k (M⁻¹)_kl α_l.
3. **Σα = 1 and boundedness.** Advect all N fractions with the same face volume fluxes that close the pressure equation, with the bounded sub-stepping of MP1 (D293 (1)), and book Σα − 1 and the clipped volume as health numbers. Solving N − 1 phases and taking the last as the remainder hides all errors in one phase, which goes negative first. (Design proposal, no source.)
4. **Closures on the total dispersed fraction.** Swarm corrections of drag, turbulent dispersion and bubble-induced turbulence must use the total gas fraction, not that of one velocity group (Frank et al. 2005, Sec. 2.4). With N phases the force registry takes pairs (k, l) and may read the sum of all dispersed phases.

Cost: each extra phase adds four transport equations (α, three momentum components) plus its scalars; the pressure equation stays one. Frank et al. (2005, Sec. 2) report that fully inhomogeneous models with one phase per size group need 15 to 20 groups and are "very expensive", which is why sizes go into groups inside a phase (section 5.3).

### 4.3 Resolved and dispersed interfaces together

A boiler drum or a separator has a free surface (resolved, VOF) and bubbles below it (dispersed). The morphology-adaptive multifield model treats the continuous and dispersed forms of the same material as separate Euler fields, resolves large interfaces "in a volume-of-fluid-like manner" through an interfacial drag, and moves mass between morphologies (HZDR project page; Schlegel et al. 2022). In Foehncraft this joins MP1 and MP3 and needs N ≥ 3 (water continuous, water dispersed drops or gas continuous, gas dispersed bubbles), which is one more reason for an N-phase core.

## 5. Several bubble sizes, also per phase

### 5.1 Why one size is not enough

The Tomiyama lift coefficient changes sign with the deformed bubble size: about 5.8 mm for air in water (Krepper et al. 2007, Sec. 3.1, Fig. 1), 5.6 mm in Li et al. (2019, Sec. 2.2). Small bubbles go to the wall, large ones to the core, and a model with one size or one gas velocity cannot produce the double peak (Krepper et al. 2007, Sec. 4.1; Li et al. 2019, Sec. 3.2, Fig. 7). In boiling the sizes also change by evaporation and condensation, so bubbles can shrink to zero size (a loss through the lower size boundary). The critical diameter depends on σ and Δρ through the Eötvös number, so for steam at high pressure it is not 5.8 mm.

### 5.2 Options per dispersed phase

| Size model | Extra fields per phase | Size-dependent velocity | Phase change | Maturity | Sources |
|---|---|---|---|---|---|
| One size (constant or correlation) | none | no | correlation only | MP3 today | |
| Interfacial area transport, one or two groups | 1 or 2 scalars | two groups only | yes (source terms) | nuclear thermal hydraulics | Ishii and co-workers, one- and two-group IATE (Purdue publication pages found, papers not read) |
| Size groups, homogeneous MUSIG | M fractions f_i | no (one velocity) | as drift between groups | industrial standard | Frank et al. 2005 Sec. 1 |
| Size groups, inhomogeneous MUSIG | M_j fractions in each of N velocity groups | yes, by velocity group | as drift | industrial standard, HZDR validation | Frank et al. 2005 Eqs. 8 to 12; Krepper et al. 2007 Eqs. 8 and 9 |
| Moments, QMOM / EQMOM / GQMOM | 2N (+1) moments | no | EQMOM handles loss at zero size | research to industrial | McGraw 1997; Yuan, Laurent, Fox 2012; Fox, Laurent, Passalacqua 2023 (via QMOMKit Sec. 1) |
| Moments, DQMOM | N weights + N abscissas | yes, one velocity per node | yes | research | Marchisio and Fox 2005 |
| Conditional moments, CQMOM / CHyQMOM | joint size-velocity moments | yes, size-conditioned | yes | research | Yuan and Fox 2011; Heylmun, Fox, Passalacqua 2021 |
| Parcels (Euler-Lagrange) | none on the mesh | yes, per parcel | yes | our L12 | de Wit et al. 2026; Radhakrishnan and Bryngelson 2026 compare EE and EL |

### 5.3 Size groups (class method), the inhomogeneous MUSIG form

Frank, Zwart, Shi, Krepper, Lucas, Rohde, "Inhomogeneous MUSIG model: a population balance approach for polydispersed bubbly flows", Int. Conf. Nuclear Energy for New Europe, Bled, 2005 (read, OSTI copy), Sec. 2.3:

- velocity group α: ∂(ρ_d r_α)/∂t + ∂(ρ_d r_α U_α,j)/∂x_j = S_α (Eq. 8)
- size group i within it: ∂(ρ_d r_α f_α,i)/∂t + ∂(ρ_d r_α f_α,i U_α,j)/∂x_j = S_α,i (Eq. 9)
- r_d = Σ_α r_α = Σ_i r_i, Σ f_i = 1 within each group, r_l + r_d = 1 (Eq. 10)
- S_α,i = B_i,B − D_i,B + B_i,C − D_i,C with Luo and Svendsen breakup and Prince and Blanch coalescence (Eq. 11); S_α = Σ_i S_α,i and Σ_α S_α = 0 (Eq. 12): breakup and coalescence move mass between groups, also across velocity groups.

Krepper, Frank, Lucas, Prasser, Zwart, ICMF 2007, Leipzig, paper 375 (read, core.ac.uk copy) give the same form (Eqs. 5 to 9) and the validation: MT-Loop test 118 with 2 velocity groups and 12 plus 22 size groups, TOPFLOW with 2 groups and 20 size classes up to 60 mm (Sec. 5.2); the breakup and coalescence rates needed the factors F_B = 0.25 and F_C = 0.05 (Sec. 5.2.1); three or four velocity groups gave no better agreement and were harder to converge (Sec. 5.1, 5.2.1); the authors call breakup and coalescence the weak points (Sec. 6). Frank et al. 2005 used 3 groups of 7 classes on TOPFLOW-074 (Sec. 3.2).

Two lessons: two velocity groups split at the lift reversal are enough for pipes; and the kernels, not the numerics, limit the accuracy. That is where the DNS of section 3 comes in.

Lehnigk et al. (AIChE J. 2022) generalise this to any number of velocity groups and add phase change as a drift between classes; I could not read the paper (Wiley and HZDR refused the connection), so I do not quote its equations.

### 5.4 Moments

From the QMOMKit survey (Denner et al. 2026, Sec. 1, read) with its references: QMOM (McGraw, Aerosol Sci. Technol. 27, 1997) takes the first 2N moments to an N-node quadrature; DQMOM (Marchisio and Fox, J. Aerosol Sci. 36, 2005) transports weights and abscissas, but "finite-volume schemes applied to the weights and abscissas do not guarantee the conservation of the moments" (citing Li, Li, Gao, Chin. J. Chem. Eng. 27, 2019); EQMOM (Yuan, Laurent, Fox, J. Aerosol Sci. 51, 2012) uses smooth kernels and so gives the population loss through a phase-space boundary "as particles or droplets shrink to zero size", which is our condensing bubble; GQMOM (Fox, Laurent, Passalacqua, J. Aerosol Sci. 167, 2023) closes recurrence coefficients and gives more nodes from the same moments; CQMOM (Yuan and Fox, JCP 230, 2011) handles joint size and velocity. Li, Marchisio, Hasse, Lucas, "Comparison of Eulerian QBMM and classical Eulerian-Eulerian method for the simulation of polydisperse bubbly flows", AIChE J. 65 (2019) e16732, CC BY (read up to Sec. 3.2): with a velocity polynomial in size (Eq. 4) the moment method reproduces the double void peak that one-size Euler-Euler cannot (Fig. 7, from the summary part), but it is "computationally demanding" and needs a "realizable finite-volume scheme" (Table 2).

Moments are the more elegant and the more fragile route (realizability, moment inversion per cell, no direct size histogram to check against). Size groups are robust, easy to check (a histogram per cell) and map one to one onto the measured distributions.

## 6. Recommendation and fit with MP3

**A phase is (material, morphology, velocity field).** Morphology: continuous, dispersed, or resolved (VOF). Two phases may share a material (water as liquid film and as drops; air as continuous gas and as bubbles of two velocity groups). Every dispersed phase carries a size model from section 5.2 behind one interface: per cell the size classes with their diameters and fractions, from which forces, heat transfer and interfacial area are evaluated per class.

**Implementation path (decided in D373 to D376):**

| Stage | Content | Checks | Size |
|---|---|---|---|
| MP3, MP4 | as decided (two phases, one size; wall boiling) | as planned | in progress |
| MP5 N-phase core | phase arrays on host and device, no gas/liquid names in the code; pairwise force registry; N×N partial elimination per cell; N-phase volume-continuity pressure; all N α advected and Σα − 1 booked; closures on the total dispersed fraction; case `[[phase]]` with any count | N = 2 bitwise against MP3; three-phase checks (two bubble velocity groups split at the lift reversal; gas, oil and water layered settling) | code about 1 to 2 hours, runs about 1 to 2 hours, CI separate |
| MP6 size groups | class method per dispersed phase (fixed classes, mass and number consistent), breakup and coalescence kernels as a registry (Luo and Svendsen, Prince and Blanch first; Morón et al. 2026 as the DNS-based option), mass exchange across velocity groups (Eq. 12 of Frank et al.), phase change as drift between classes and loss below the smallest class | homogeneous breakup and coalescence against analytic solutions; MT-Loop-type pipe with two velocity groups; size-distribution slopes around d_H as a sanity check | code about 2 to 3 hours, runs about 1 to 2 hours |
| MP7 moments | EQMOM or GQMOM as an alternative size model behind the same interface; CQMOM for size-conditioned velocity | the MP6 cases, moment against class results | not sized yet |
| MP8 resolved and dispersed together | morphology-adaptive coupling of MP1 VOF and the N-phase core | rising bubble, drum free surface with bubbles | not sized yet |
| later | interface-resolved N-phase DNS (conservative phase field after Huang et al.), graphics benchmark of section 4.1 | | not sized yet |

**Agents first (D77 to D87).** Every phase and every size class appears in the structured output: per-phase volume and mass ledgers, Σα − 1, size histograms per region, d32 fields, breakup and coalescence rates, so an agent can check a polydisperse run without a picture.

## 7. Requirement text (added to L3)

L3 amendment: "Multiphase with any number of phases. A phase is a velocity field of one material in one morphology (continuous, dispersed, resolved); the code has no fixed limit on the number of phases (arrays per phase, pairwise closures), practical cases use two to five. Each dispersed phase has its own size model: one size, size groups (class method; several velocity groups of one material exchange mass by breakup and coalescence), or moments; phase change moves mass between size classes. The two-phase rows of MP3 stay bit for bit as the case N = 2." The formulation of each stage gets its D-entries when the stage is built, as D293 and D351 did for MP1 and MP2.

## 8. Decisions (answered 2026-11-28)

| ID | Question | Options | Answer |
|---|---|---|---|
| D373 | When does the multiphase core become general for N phases? | (a) MP5, right after MP3 and MP4, two-phase rows bit for bit as N = 2; (b) MP3 builds N phases now; (c) only when a case needs it | **(a)** |
| D374 | Which bubble size model do dispersed phases get first? | (a) size groups first (MP6), moments later (MP7); (b) moments first (EQMOM or GQMOM); (c) interfacial area only (two groups); (d) bubbles as parcels (L12) | **(a)** |
| D375 | How do bubbles of different sizes get different velocities? | (a) velocity groups as phases of one material (inhomogeneous MUSIG); (b) size-conditioned velocity by conditional moments; (c) one velocity per phase with an algebraic slip per class | **(a)** |
| D376 | When do free surfaces and dispersed bubbles share one model? | (a) morphology-adaptive hybrid as MP8 after the size groups, interface-resolved N-phase DNS a later item; (b) hybrid right after MP5; (c) no hybrid, DNS later | **(a)** |

The stage plan of section 6 is therefore the decided order: MP5 N-phase core, MP6 size groups with velocity groups as phases, MP7 moments, MP8 resolved and dispersed interfaces together, then the interface-resolved N-phase DNS with the graphics benchmark of section 4.1.

## 9. Sources

Read (full or in large part): arXiv:2406.01322; arXiv:2309.13995; arXiv:2309.13976; CPC 316 (2025) 109804 (first half); arXiv:2010.01099 (J. Comput. Appl. Math. 408, 2022); arXiv:2101.04252 (part); Frank et al. NENE 2005 (OSTI); Krepper et al. ICMF 2007 (core.ac.uk); Li et al. AIChE J. 65 (2019) e16732 (to Sec. 3.2); arXiv:2610.36189 (Sec. 1 and references); abstracts and first pages of arXiv:2605.07504, 2402.19121, 2510.03968, 2607.27426, 2511.21031, 2603.13697, 2112.14172; Keenan Crane's publication list; Ren et al. TOG 2014.

## 10. What I could not read

- Soligo, Roccon, Soldati, J. Fluids Eng. 2021 and J. Fluid Mech. 881 (2019); Roccon, Zonta, Soldati, Phys. Rev. Fluids 2023 (only abstract and metadata).
- Lehnigk et al., AIChE J. 2022 (Wiley and HZDR refused).
- Yuan, Laurent, Fox 2012 (EQMOM) and the GQMOM paper (HAL bot check); they are cited here through QMOMKit and Li et al.
- Karema and Lo 1999, Spalding's partial elimination, Darwish and Moukalled's multi-fluid papers (not open); Hanimann et al. 2021 only title, abstract and a summary.
- Schlegel et al. Nucl. Sci. Eng. 2022 (only the project page).
- The one- and two-group interfacial area transport papers of Ishii and co-workers (only publication pages).
- Kumar and Ramkrishna 1996 (fixed pivot), Luo and Svendsen 1996, Prince and Blanch 1990 (cited through Frank et al. and Krepper et al.).
- Krepper and Rzehak, subcooled boiling with MUSIG (HZDR pages refused).
- The MP3 code (PR not open yet): section 4.2 is to be checked against it when it is.
