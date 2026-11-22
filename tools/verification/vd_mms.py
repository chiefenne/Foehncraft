"""The manufactured solution of the variable density (stage 7c, check V7c.4).

Writes the momentum and heat sources of `[analytic] type =
"variable-density"` (crates/foehncraft-flow/src/analytic.rs) into the case
examples/cases/flow/vd_mms_steady.toml. On the unit square, an ideal gas
with R/M = 1 at p₀ = 1 (so ρ = 1/T), c_p = 3.5 (γ = 1.4), constant μ and k:

    T   = 1 + a sin πx sin πy
    ψ   = A sin²πx sin²πy,  ρu = ∂ψ/∂y,  ρv = −∂ψ/∂x  (∇·(ρu) = 0 exactly)
    p   = P cos πx cos πy

    f   = ∇·(ρ u u) + ∇p − ∇·τ,  τ = μ (∇u + ∇uᵀ − ⅔ (∇·u) I)
    q   = c_p ρu·∇T − k ∇²T                      (steady: no pressure work)

    pip install sympy
    python3 -m tools.verification.vd_mms            write the case
    python3 -m tools.verification.vd_mms --check    fail when it is out of date
"""
import pathlib
import sys

import sympy as sp

A, AT, P = 0.02, 0.5, 0.01
MU, K, CP = 0.002, 0.007, 3.5

CASE = pathlib.Path(__file__).resolve().parents[2] / "examples/cases/flow/vd_mms_steady.toml"


def sources():
    x, y = sp.symbols("x y")
    pi = sp.pi
    t = 1 + AT * sp.sin(pi * x) * sp.sin(pi * y)
    rho = 1 / t
    psi = A * sp.sin(pi * x) ** 2 * sp.sin(pi * y) ** 2
    mx, my = sp.diff(psi, y), -sp.diff(psi, x)
    u, v = mx / rho, my / rho
    p = P * sp.cos(pi * x) * sp.cos(pi * y)
    div = sp.diff(u, x) + sp.diff(v, y)
    g = [[sp.diff(u, x), sp.diff(u, y)], [sp.diff(v, x), sp.diff(v, y)]]
    tau = [
        [MU * (g[i][j] + g[j][i] - (sp.Rational(2, 3) * div if i == j else 0)) for j in range(2)]
        for i in range(2)
    ]
    vel = [u, v]
    flux = [mx, my]
    f = []
    for i in range(2):
        conv = sp.diff(flux[0] * vel[i], x) + sp.diff(flux[1] * vel[i], y)
        visc = sp.diff(tau[i][0], x) + sp.diff(tau[i][1], y)
        f.append(conv + sp.diff(p, [x, y][i]) - visc)
    q = CP * (mx * sp.diff(t, x) + my * sp.diff(t, y)) - K * (sp.diff(t, x, 2) + sp.diff(t, y, 2))
    return [text(e) for e in f], text(q)


def text(e):
    """An expression in the syntax of the case files."""
    return sp.sstr(e).replace("**", "^")


def case():
    (fx, fy), q = sources()
    return f'''# The manufactured solution of the variable density (stage 7c, D151,
# check V7c.4): written by `python3 -m tools.verification.vd_mms`; do not
# edit the source expressions by hand.
#
#   cargo run --release -p foehncraft-cli -- flow examples/cases/flow/vd_mms_steady.toml
#
# On the unit square, an ideal gas with R/M = 1 (molar mass 8.314 kg/mol)
# at the operating pressure p₀ = 1 Pa, so ρ = 1/T; c_p = {CP} J/(kg K)
# (γ = 1.4), μ = {MU} Pa s and k = {K} W/(m K) constant. The exact solution
# (`[analytic] type = "variable-density"`): T = 1 + a sin πx sin πy with
# a = {AT}, the mass flux of the stream function ψ = A sin²πx sin²πy with A
# = {A} (∇·(ρu) = 0 exactly, u = 0 on the walls), and p = P cos πx cos πy
# with P = {P}. The momentum source (`[flow] body_force`) and the heat
# source (`[energy] source`) make it a solution of the steady low-Mach
# equations with the full viscous stress μ (∇u + ∇uᵀ − ⅔ (∇·u) I). The
# box is closed, so p₀ follows from the mass of the start temperature
# (the exact one). Three nested meshes (16², 32², 64²) with the coupled
# algorithm (SIMPLEC gives the same solution in more iterations): the L2
# errors of u and T fall at second order (observed 1.92 and 1.90 on the
# finest pair; on the polygonal mesh of `type = "poly"` 2.32 and 2.16).
# On the randomly perturbed hexahedra of `perturb` the velocity order
# drops on the finest pair with a constant density as well (a property of
# the perturbed meshes, not of the variable density).

case_format = "0.6"

[mesh.generate]
type = "hex"
dimension = 2
upper = [1.0, 1.0, 0.05]
cells = [16, 16, 1]

[materials.gas]
kind = "gas"
description = "an ideal gas with R/M = 1 J/(kg K) (manufactured solution)"
molar_mass = 8.314462618
density = {{ form = "ideal-gas", source = "ideal gas, rho = p M / (R T)" }}
specific_heat = {CP}
conductivity = {K}
viscosity = {MU}

[fluid]
material = "gas"
temperature = 1.0
operating_pressure = 1.0

[flow]
body_force = [
  "{fx}",
  "{fy}",
  0.0,
]

[energy]
initial = "1 + {AT} * sin(pi * x) * sin(pi * y)"
exact = "1 + {AT} * sin(pi * x) * sin(pi * y)"
source = "{q}"
density = "variable"
thermodynamic_pressure = "closed-volume"

[boundary.xmin]
type = "wall"
[boundary.xmin.energy]
temperature = 1.0
[boundary.xmax]
type = "wall"
[boundary.xmax.energy]
temperature = 1.0
[boundary.ymin]
type = "wall"
[boundary.ymin.energy]
temperature = 1.0
[boundary.ymax]
type = "wall"
[boundary.ymax.energy]
temperature = 1.0
[boundary.zmin]
type = "symmetry"
[boundary.zmax]
type = "symmetry"

[solver]
algorithm = "coupled"
tolerance = 1e-10
max_iterations = 3000

[solver.energy]
tolerance = 1e-12
max_iter = 200

[analytic]
type = "variable-density"
flux = {A}
temperature = {AT}
pressure = {P}

[verification]
levels = [0, 1, 2]
expected_order = [1.8, 2.2]
scalar_order = [1.8, 2.2]
require_convergence = true

[output]
vtkhdf = false
wall_csv = false
'''


def main():
    text_ = case()
    if "--check" in sys.argv[1:]:
        if CASE.read_text() != text_:
            sys.exit(f"{CASE} is out of date: run python3 -m tools.verification.vd_mms")
        return
    CASE.write_text(text_)
    print(f"wrote {CASE}")


if __name__ == "__main__":
    main()
