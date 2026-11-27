"""Similarity solution of laminar forced convection along an isothermal
flat plate (Pohlhausen 1921), the reference of the wall heat quantities of
`foehn post derive` (case plate_heated_laminar).

With η = y √(U∞/(ν x)), the stream function ψ = √(ν U∞ x) f(η) of Blasius
(1908) and θ = (T − T∞)/(T_w − T∞) = θ(η):

    f''' + ½ f f''  = 0,   f(0) = f'(0) = 0,  f'(∞) = 1,
    θ''  + ½ Pr f θ' = 0,  θ(0) = 1,          θ(∞) = 0.

The energy equation integrates once: θ'(η) = θ'(0) exp(−½ Pr ∫₀^η f), so
−θ'(0) = 1 / ∫₀^∞ exp(−½ Pr ∫₀^η f dη') dη, and the local Nusselt number
is Nu_x = h x / k = −θ'(0) Re_x^(1/2). Pohlhausen's fit of it is the
correlation Nu_x = 0.332 Re_x^(1/2) Pr^(1/3) (0.6 < Pr < 10). The Blasius
wall value f''(0) = 0.332057 comes from shooting (bisection on the far
field); the script checks it.

    python tools/verification/pohlhausen.py             # Pr 0.71
    python tools/verification/pohlhausen.py --pr 0.71 --profile out.csv
"""

from __future__ import annotations

import argparse

import numpy as np


def _rk4(s2: float, eta_max: float, n: int, pr: float):
    """Integrates Blasius and the integral ∫ f from η = 0 with f''(0) = s2:
    η, f, f', f'', ∫₀^η f."""
    h = eta_max / n

    def rhs(s):
        f, f1, f2, g = s
        return np.array([f1, f2, -0.5 * f * f2, f])

    s = np.array([0.0, 0.0, s2, 0.0])
    out = np.empty((n + 1, 5))
    out[0] = [0.0, *s]
    for i in range(n):
        k1 = rhs(s)
        k2 = rhs(s + 0.5 * h * k1)
        k3 = rhs(s + 0.5 * h * k2)
        k4 = rhs(s + h * k3)
        s = s + h / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        out[i + 1] = [(i + 1) * h, *s]
    return out


def blasius_wall(eta_max: float = 12.0, n: int = 12000) -> float:
    """f''(0) of Blasius, f'(η_max) = 1 by bisection."""
    lo, hi = 0.1, 1.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if _rk4(mid, eta_max, n, 0.0)[-1, 2] > 1.0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def solve(pr: float, eta_max: float = 12.0, n: int = 12000):
    """(f''(0), −θ'(0), table η f f' θ) at Prandtl number `pr`."""
    s2 = blasius_wall(eta_max, n)
    t = _rk4(s2, eta_max, n, pr)
    eta, g = t[:, 0], t[:, 4]
    w = np.exp(-0.5 * pr * g)
    # trapezoidal integrals of exp(−½ Pr ∫ f): the denominator and θ
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (w[1:] + w[:-1]) * np.diff(eta))])
    total = cum[-1]
    theta = 1.0 - cum / total
    return s2, 1.0 / total, np.column_stack([eta, t[:, 1], t[:, 2], theta])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--pr", type=float, default=0.71)
    ap.add_argument("--profile", help="write η, f, f', θ to this CSV")
    a = ap.parse_args()
    s2, q, table = solve(a.pr)
    assert abs(s2 - 0.332057) < 2e-6, s2
    print(f"Blasius f''(0) = {s2:.6f}")
    print(f"Pr = {a.pr}: -theta'(0) = {q:.6f}, Nu_x / Re_x^(1/2) = {q:.6f} "
          f"(correlation 0.332 Pr^(1/3) = {0.332 * a.pr ** (1 / 3):.6f})")
    if a.profile:
        np.savetxt(a.profile, table, delimiter=",", header="eta,f,f1,theta", comments="")


if __name__ == "__main__":
    main()
