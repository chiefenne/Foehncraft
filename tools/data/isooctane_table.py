"""Write a saturation table of liquid iso-octane (2,2,4-trimethylpentane) in
the column layout of the NIST WebBook tables, for liquid_properties.py (D320).

Iso-octane is not in the NIST fluid system, so the table comes from published
correlations (only our fit of them goes into the repository):

    vapour pressure  Wagner 3-6 form, McGarry, Ind. Eng. Chem. Process Des. Dev.
                     22 (1983) 313 (Tc 543.96 K, pc 2.56155 MPa, from 265 K);
                     checked against Willingham, Taylor, Pignocco and Rossini,
                     J. Res. NBS 35 (1945) 219 (Antoine, 298-373 K)
    density          DIPPR 105, Perry's Chemical Engineers' Handbook 8th ed.,
                     Table 2-32
    specific heat    DIPPR 100, Perry 8th ed., Table 2-153 (to 520 K)
    latent heat      DIPPR 106, Perry 8th ed., Table 2-150; checked against the
                     calorimetric values of Majer and Svoboda (NIST WebBook)
    viscosity        DIPPR 101, Perry 8th ed., Table 2-313
    surface tension  Mulero, Cachadina and Parra, J. Phys. Chem. Ref. Data 41
                     (2012) 043105 (273-373 K, extrapolated above)

Above about 400 K the density, specific heat, viscosity and surface tension
rest on the correlations' extrapolation (a few per cent, more near Tc).

    python3 tools/data/isooctane_table.py DIR     # writes DIR/iso-octane.txt
"""

import math
import sys

M = 0.11423  # kg/mol
TC_W, PC_W = 543.96, 2.56155e6
TC = 543.8  # the DIPPR end, the last row of the table


def psat(t):
    tr = t / TC_W
    s = 1.0 - tr
    return PC_W * math.exp((-7.3889 * s + 1.25294 * s**1.5 - 3.16606 * s**3 - 2.22001 * s**6) / tr)


def density(t):
    return 0.59059 / 0.27424 ** (1.0 + (1.0 - t / TC) ** 0.2847) * M * 1e3


def cp(t):
    t = min(t, 520.0)
    return (95275 + 696.7 * t - 1.3765 * t * t + 0.0021734 * t**3) / (M * 1e3)


def latent(t):
    return 47.711e6 * (1.0 - t / TC) ** 0.37949 / (M * 1e3)


def viscosity(t):
    return math.exp(-12.928 + 1137.5 / t + 0.25725 * math.log(t) - 3.6929e-28 * t**10)


def sigma(t):
    return 0.04794 * max(1.0 - t / 543.87, 0.0) ** 1.209


def main(argv):
    head = [
        "Temperature (K)", "Pressure (MPa)", "Density (l, kg/m3)", "Enthalpy (l, kJ/kg)",
        "Cp (l, J/g*K)", "Viscosity (l, Pa*s)", "Surf. Tension (l, N/m)",
        "Density (v, kg/m3)", "Enthalpy (v, kJ/kg)",
    ]
    rows = []
    t = 265.0
    while t < TC - 0.5:
        rows.append([t, psat(t) * 1e-6, density(t), 0.0, cp(t) * 1e-3, viscosity(t),
                     sigma(t), 0.0, latent(t) * 1e-3])
        t += 1.0
    rc = 0.59059 / 0.27424 * M * 1e3
    rows.append([TC, psat(TC) * 1e-6, rc, 0.0, float("nan"), float("nan"), 0.0, rc, 0.0])
    with open(f"{argv[1]}/iso-octane.txt", "w") as f:
        f.write("\t".join(head) + "\n")
        for r in rows:
            f.write("\t".join(f"{x:.6g}" for x in r) + "\n")


if __name__ == "__main__":
    main(sys.argv)
