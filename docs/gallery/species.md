# Species and reacting flows

Part of the [Example gallery](../gallery.md): mixtures of ideal gases whose species are carried by the flow, the first stage (SF1) of the [development path to an industrial steam furnace](../assessments/steam-furnace-development-path.md). Chemistry and heat release follow in SF2 and SF3.

| # | Example | Stage | Feature |
|---|---|---|---|
| 278 | [A methane jet in a coflow of air, and binary diffusion against the exact solution](#species-transport) | SF1 | `[species]`, `names`, `closure`, `initial`, `diffusion`, `schmidt`, `lewis`, `turbulent_schmidt`, `properties`, `[boundary.<name>.species]`, fields `Y_<species>` |

---

<a id="species-transport"></a>
## 278. A methane jet in a coflow of air, and binary diffusion against the exact solution

<img src="../figures/gallery/species_coaxial_jet.png" alt="Three panels of an axisymmetric methane jet 300 mm long, mirrored about its axis: the methane mass fraction falls from 1 at the nozzle of radius 5 mm to about 0.55 at the far end as the jet widens; the density of the mixture is 0.65 kg/m³ in the jet core and 1.17 kg/m³ in the air around it; the speed is 1 m/s in the jet core and 0.2 m/s in the coflow. Bottom: along the axis the methane mass fraction stays at 1 for the first 50 mm and then falls, oxygen and nitrogen rise, and the density rises from 0.65 to 0.82 kg/m³." width="1200" height="960">

`species_coaxial_jet.toml` blows pure methane from a nozzle of radius
5 mm at 1 m/s into a coflow of air at 0.2 m/s, both at 300 K and
1 atm (Re = 590, laminar), on an axisymmetric wedge of 2°. The mixture has
three species: CH₄ and O₂ are transported, N₂ is the closure species,
Y_N2 = 1 − Y_CH4 − Y_O2 in every cell. The jet is isothermal, yet its
density varies by a factor of 1.8: ρ = p₀ M/(R T) with the molar mass of
the local mixture, 1/M = Σ Y_i/M_i, through the low-Mach variable density
of stage 7c. c_p, λ and μ are those of the local mixture: NASA-7
polynomials and kinetic theory (Chapman–Enskog with the collision
integrals of Neufeld et al. 1972, the Wilke and the Mathur mixing rules),
with the 53 species of GRI-Mech 3.0. Diffusion: ρ D_i = μ/Sc, Sc = 0.7.
The run converges in 154 outer iterations (9 s on one core). The carbon,
hydrogen, oxygen and nitrogen accounts of the ledger close to 6.4·10⁻⁹ of
their inflow (the solver tolerance 10⁻⁸).

<img src="../figures/gallery/species_binary_diffusion.png" alt="Mass fraction of methane in a slab 1 m long: the start cosine from 0.9 to 0.1 (dotted), the exact solution after 50 s (black line) from 0.74 to 0.26, the computed values on it (green circles), and nitrogen as one minus methane (grey squares)." width="800" height="440">

`species_binary_diffusion.toml` checks the transport against an exact
solution: CH₄ in N₂ at rest in a closed slab, Y_CH4 = 0.5 + 0.4 cos(πx) at
the start, D = 10⁻³ m²/s. The mode decays as exp(−Dπ²t); after 50 s with
BDF2 on 100 cells the error is 9.5·10⁻⁶ (L2), 1.3·10⁻⁵ at worst. The
closure fraction is 1 − Y_CH4 exactly, and the accounts of both species
and of C, H and N close to 10⁻¹⁹ kg/s.

**How it is built.** Each transported species is a scalar field `Y_<name>`
of the generic transport equation (D116), solved with the density of the
flow in its time term and its convection, so the device kernels of the
scalars (CPU and GPU, one source) and their path on several ranks carry
it. The convection is bounded by the vertex limiter (MLP, D237); the mass
fractions are never clipped after a solve. Once per outer iteration the
mixture is evaluated per cell and face from the new mass fractions and
temperature (`foehncraft_kernels::species`, a cell-local body with a
fixed summation order): it gives the face diffusivities of the species,
c_p, λ and μ of the energy equation and the flow, and the molar mass of
the density. The ledger (D119) gets an account per species and per
chemical element: convected and diffused mass per boundary, and the
continuity defect Σ Y_P (∇·ṁ)_P V_P that the non-conservative form of the
equation carries, so each account closes to the solver tolerance.

**What this version does not do.** No chemistry and no heat release
(SF2, SF3), no enthalpy flux of the species diffusion Σ h_i J_i in the
energy equation (it matters with heat release; zero in an isothermal
mixture), no mixture-averaged or multicomponent diffusion with a
correction velocity (differential diffusion of hydrogen, SF2h): every
species has a Schmidt or a Lewis number. One fluid region: no solids,
porous regions or several fluids with species (SPE002). On ranks the variable density (7c) and the MLP limiter are
refused, so a species case runs there with a constant density. The mixture is
evaluated on the host from the device fields once per outer iteration
(the same round trip as the property refresh of the energy equation); a
device kernel of it follows with the chemistry (SF2). The MLP limiter
lets smooth extrema of ε = 1 % of the range pass, so a mass fraction may
leave [0, 1] by a few 10⁻⁶ (the jet: −5.9·10⁻⁷ and 1 + 3·10⁻⁶); a sum
above 1 + 10⁻⁴ is warned about (SPE003). The ledger takes the
linear-upwind reconstruction without the clip of `bounded = true`.

**What is checked.** The exact decay of binary diffusion (L2 error below
10⁻⁴); the element accounts of the jet (below 10⁻⁶ of the inflow) and the
jet density p₀ M_CH4/(R T); pure-gas viscosities of N₂, O₂, CH₄, Ar, H₂
and CO₂ at 300 K within 3 % of measured values and conductivities within
5 %; air as a mixture (M, c_p, μ, λ); the enthalpy of formation of CH₄; the
binary diffusion coefficient of CH₄–N₂; the species case (constant density, mixture
properties) on 1 to 4 ranks bitwise equal to one rank; the refusals SPE001.

Files: `fields.vtkhdf` (`Y_CH4`, `Y_O2`, `density`), `line_axis.csv`,
`line_x50mm.csv`, `result.json` (`balance.last.ledger`: `species …`,
`element …`). [Recipe 278](../recipes/species.md#species-transport),
cases in the [catalogue](../cases/flow.md#species), images by
`python3 -m tools.plot.species_fields`.
