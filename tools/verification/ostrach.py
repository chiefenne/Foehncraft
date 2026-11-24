"""Similarity solution of laminar free convection along an isothermal
vertical plate (Ostrach 1953), the reference of case V7.12 (stage 7b2).

With η = (y/x) (Gr_x/4)^(1/4), the stream function
ψ = 4 ν (Gr_x/4)^(1/4) F(η) and θ = (T − T∞)/(T_w − T∞) = H(η):

    F''' + 3 F F'' − 2 F'² + H = 0,
    H''  + 3 Pr F H'          = 0,
    F(0) = F'(0) = 0, H(0) = 1,   F'(∞) = H(∞) = 0.

The local Nusselt number is Nu_x = −H'(0) (Gr_x/4)^(1/4) and the velocity
along the plate u = (2 ν/x) Gr_x^(1/2) F'(η). The two unknown wall values
F''(0), H'(0) are found by shooting (Newton on the far-field conditions,
the far edge moved out until the values settle). Ostrach's table gives
−H'(0) = 0.5046 at Pr = 0.72; the script checks that value.

    python tools/verification/ostrach.py            # Pr 0.71 and 0.72
    python tools/verification/ostrach.py --pr 0.71 --profile out.csv
"""

from __future__ import annotations

import argparse

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve


def _rhs(pr: float):
    def f(_eta, s):
        F, F1, F2, H, H1 = s
        return [F1, F2, -3.0 * F * F2 + 2.0 * F1 * F1 - H, H1, -3.0 * pr * F * H1]

    return f


def _shoot(pr: float, guess, eta_max: float):
    def miss(g):
        sol = solve_ivp(
            _rhs(pr), (0.0, eta_max), [0.0, 0.0, g[0], 1.0, g[1]], rtol=1e-12, atol=1e-14
        )
        end = sol.y[:, -1]
        return [end[1], end[3]]

    g, _info, ok, msg = fsolve(miss, guess, xtol=1e-13, full_output=True)
    if ok != 1:
        raise RuntimeError(f"shooting did not converge at eta_max {eta_max}: {msg}")
    return g


def solve(pr: float):
    """Returns (F''(0), H'(0)) for Prandtl number pr."""
    g = np.array([0.68, -0.5])
    for eta_max in (6.0, 8.0, 10.0, 12.0, 14.0):
        g = _shoot(pr, g, eta_max)
    return float(g[0]), float(g[1])


def profile(pr: float, eta_max: float = 10.0, n: int = 201):
    """η, F'(η), H(η) on n points."""
    f2, h1 = solve(pr)
    eta = np.linspace(0.0, eta_max, n)
    sol = solve_ivp(
        _rhs(pr), (0.0, eta_max), [0.0, 0.0, f2, 1.0, h1], t_eval=eta, rtol=1e-12, atol=1e-14
    )
    return eta, sol.y[1], sol.y[3]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--pr", type=float, action="append")
    ap.add_argument("--profile", help="write η, F', H to this CSV file")
    a = ap.parse_args()
    for pr in a.pr or [0.71, 0.72]:
        f2, h1 = solve(pr)
        print(f"Pr {pr}: F''(0) = {f2:.6f}, -H'(0) = {-h1:.6f}")
    if a.pr is None:
        _f2, h1 = solve(0.72)
        assert abs(-h1 - 0.5046) < 5e-4, "Ostrach's table: -H'(0) = 0.5046 at Pr 0.72"
        print("check: -H'(0) at Pr 0.72 matches Ostrach's 0.5046")
    if a.profile:
        eta, fp, h = profile((a.pr or [0.71])[0])
        np.savetxt(a.profile, np.c_[eta, fp, h], delimiter=",", header="eta,Fp,H", comments="")


if __name__ == "__main__":
    main()
