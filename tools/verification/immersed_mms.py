"""The manufactured solution of the penalised equations (L13 stage 1, D226).

Writes the momentum source of `[analytic] type = "variable-density"` with a
zero temperature amplitude (so a constant density ρ = 1) and a smooth solid
fraction into the case examples/cases/flow/immersed_mms.toml. On the unit
square, μ constant, the permeability K:

    χ   = exp(−((x − ½)² + (y − ½)²) / w)      the solid fraction (1 at the middle)
    ψ   = A sin²πx sin²πy,  u = ∂ψ/∂y,  v = −∂ψ/∂x   (∇·u = 0, u = 0 on the walls)
    p   = P cos πx cos πy

    f   = ρ (u·∇) u + ∇p − μ ∇²u + (μ χ / K) u

The Brinkman sink acts everywhere with a smooth χ, so the check measures the
order of the discretisation of the penalised equations (the sink inside the
momentum rows, the Rhie–Chow flux and the SIMPLEC factor that see it), not
the smeared wall of a body.

    pip install sympy
    python3 -m tools.verification.immersed_mms            write the case
    python3 -m tools.verification.immersed_mms --check    fail when it is out of date
"""
import pathlib
import sys

import sympy as sp

A, P = 0.02, 0.01
MU, KPERM, W = 0.01, 1e-3, 0.05

CASE = pathlib.Path(__file__).resolve().parents[2] / "examples/cases/flow/immersed_mms.toml"

CHI = f"exp(-((x - 0.5)^2 + (y - 0.5)^2) / {W})"


def sources():
    x, y = sp.symbols("x y")
    pi = sp.pi
    chi = sp.exp(-((x - sp.Rational(1, 2)) ** 2 + (y - sp.Rational(1, 2)) ** 2) / W)
    psi = A * sp.sin(pi * x) ** 2 * sp.sin(pi * y) ** 2
    u, v = sp.diff(psi, y), -sp.diff(psi, x)
    p = P * sp.cos(pi * x) * sp.cos(pi * y)
    f = []
    for w, d in ((u, x), (v, y)):
        conv = u * sp.diff(w, x) + v * sp.diff(w, y)
        visc = MU * (sp.diff(w, x, 2) + sp.diff(w, y, 2))
        f.append(conv + sp.diff(p, d) - visc + MU / KPERM * chi * w)
    return [text(e) for e in f]


def text(e):
    """An expression in the syntax of the case files."""
    return sp.sstr(e).replace("**", "^")


def case():
    fx, fy = sources()
    return f'''# The manufactured solution of the penalised equations (L13 stage 1, D226):
# written by `python3 -m tools.verification.immersed_mms`; do not edit the
# source expressions by hand.
#
#   cargo run --release -p foehncraft-cli -- flow examples/cases/flow/immersed_mms.toml
#
# On the unit square, ρ = 1, μ = {MU} Pa s, K = {KPERM} m² (μ/K = {MU / KPERM:g}, the
# sink at the middle about a hundred times the viscous term). A smooth solid
# fraction χ = exp(−((x − ½)² + (y − ½)²)/{W}) (a `field` body) and the exact
# solution of `[analytic] type = "variable-density"` with a zero temperature
# amplitude (a constant density): the velocity of the stream function
# ψ = A sin²πx sin²πy, A = {A} (u = 0 on the walls), p = P cos πx cos πy,
# P = {P}. The momentum source makes it a solution of the steady
# Navier–Stokes equations with the Brinkman term μ χ u / K. Three nested
# meshes (16², 32², 64²): the L2 error of u falls at second order on
# hexahedra, triangular prisms (`type = "tet"`), polygons (`poly`) and the
# mixed mesh (`mixed`), as without penalisation (tests/immersed.rs).

case_format = "0.7"

[mesh.generate]
type = "hex"
dimension = 2
upper = [1.0, 1.0, 0.05]
cells = [16, 16, 1]

[fluid]
density = 1.0
viscosity = {MU}
pressure_values = "gauge"

[flow]
body_force = [
  "{fx}",
  "{fy}",
  0.0,
]

[immersed]
permeability = {KPERM}

[[immersed.body]]
name = "density"
shape = "field"
fraction = "{CHI}"

[boundary.xmin]
type = "wall"
[boundary.xmax]
type = "wall"
[boundary.ymin]
type = "wall"
[boundary.ymax]
type = "wall"
[boundary.zmin]
type = "symmetry"
[boundary.zmax]
type = "symmetry"

[solver]
algorithm = "coupled"
tolerance = 1e-10
max_iterations = 3000

[analytic]
type = "variable-density"
flux = {A}
temperature = 0.0
pressure = {P}

[verification]
levels = [0, 1, 2]
expected_order = [1.8, 2.3]
require_convergence = true

[output]
vtkhdf = false
wall_csv = false
'''


def main():
    text_ = case()
    if "--check" in sys.argv[1:]:
        if CASE.read_text() != text_:
            sys.exit(f"{CASE} is out of date: run python3 -m tools.verification.immersed_mms")
        return
    CASE.write_text(text_)
    print(f"wrote {CASE}")


if __name__ == "__main__":
    main()
