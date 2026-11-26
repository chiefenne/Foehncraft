"""Fit the liquid property correlations of the spray liquids (D320) to the
saturation tables of the NIST Chemistry WebBook (SRD 69, "Thermophysical
Properties of Fluid Systems", Lemmon, McLinden and Friend) and print the
coefficients as the Rust data of `crates/foehncraft-kernels/src/liquids.rs`.

Only our own fitted coefficients go into the repository; the NIST tables stay
outside (download them with the URLs printed by `--urls`). Forms, with
t = 1 - T/Tc and Tr = T/Tc:

    vapour pressure  ln(p/pc) = (a t + b t^1.5 + c t^2.5 + d t^5)/Tr   (Wagner)
    density          rho = rho_c (1 + sum a_k t^e_k), e = 0.35, 2/3, 1, 4/3, 5/3, 16/3
    specific heat    c_l = a + b/t + c t + d t^2
    latent heat      L = t^0.38 (a + b t + c t^2 + d t^3)                 (Watson)
    surface tension  sigma = t^1.26 (a + b t + c t^2 + d t^3)
    viscosity        ln mu = sum a_k Tr^e_k, e = 0, -1, -2, -3, -4, 1

Density, latent heat and surface tension are fitted up to Tr = 0.98, the
specific heat and the viscosity up to Tr = 0.95, the vapour pressure to Tc;
each fit minimises the relative error.

    python3 tools/data/liquid_properties.py --urls
    python3 tools/data/liquid_properties.py DIR      # DIR holds <name>.txt tables
"""

import csv
import math
import sys

import numpy as np

FLUIDS = {
    # name: (NIST id, molar mass kg/mol, T low, T high)
    "water": ("C7732185", 0.01801528, 273.16, 647.0),
    "n-heptane": ("C142825", 0.100202, 182.6, 540.0),
    "n-dodecane": ("C112403", 0.170335, 263.6, 658.0),
    # not in the NIST fluid system: the table of isooctane_table.py
    "iso-octane": (None, 0.11423, 265.0, 543.8),
}

DENSITY_EXP = [0.35, 2 / 3, 1.0, 4 / 3, 5 / 3, 16 / 3]
VISCOSITY_EXP = [0.0, -1.0, -2.0, -3.0, -4.0, 1.0]
R = 8.314462618

URL = (
    "https://webbook.nist.gov/cgi/fluid.cgi?Action=Data&Wide=on&ID={id}&Type=SatT"
    "&Digits=6&THigh={hi}&TLow={lo}&TInc=1&RefState=DEF&TUnit=K&PUnit=MPa"
    "&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=Pa*s&STUnit=N%2Fm"
)


def num(s):
    try:
        return float(s)
    except ValueError:
        return math.nan


def read(path):
    """Rows of the tab-delimited NIST table: liquid and vapour columns side by side."""
    rows = list(csv.reader(open(path), delimiter="\t"))
    head = rows[0]
    out = {}
    for k, h in enumerate(head):
        out[h.strip()] = np.array([num(r[k]) for r in rows[1:] if len(r) > k])
    return out


def lsq(cols, y, w=None):
    """Least squares of y on the columns, weighted by w (relative error by default)."""
    a = np.column_stack(cols)
    w = 1.0 / np.abs(y) if w is None else w
    c, *_ = np.linalg.lstsq(a * w[:, None], y * w, rcond=None)
    return c, a @ c


def fit(name, path):
    d = read(path)
    t_k = d["Temperature (K)"]
    p = d["Pressure (MPa)"] * 1e6
    keys = list(d)
    rho_l = d[[k for k in keys if k.startswith("Density (l")][0]]
    rho_v = d[[k for k in keys if k.startswith("Density (v")][0]]
    h_l = d[[k for k in keys if k.startswith("Enthalpy (l")][0]] * 1e3
    h_v = d[[k for k in keys if k.startswith("Enthalpy (v")][0]] * 1e3
    cp_l = d[[k for k in keys if k.startswith("Cp (l")][0]] * 1e3
    mu_l = d[[k for k in keys if k.startswith("Viscosity (l")][0]]
    sig = d[[k for k in keys if k.startswith("Surf. Tension")][0]]
    tc, pc, rc = t_k[-1], p[-1], 0.5 * (rho_l[-1] + rho_v[-1])
    tr = t_k / tc
    t = 1.0 - tr
    res = {"tc": tc, "pc": pc, "rho_c": rc, "t_min": t_k[0]}
    report = []

    def band(mask, y, f, what):
        """Largest relative error of fit f against the table y on mask."""
        r = np.abs(f / y - 1.0)
        tt = t_k[mask]
        e9 = np.max(r[tt <= 0.9 * tc])
        report.append(f"{what} {e9:.1e} to 0.9 Tc, {np.max(r):.1e} to {tt[-1]:.0f} K")

    m = t > 1e-4
    c, f = lsq([t[m], t[m] ** 1.5, t[m] ** 2.5, t[m] ** 5], np.log(p[m] / pc) * tr[m], np.ones(m.sum()))
    res["psat"] = c
    band(m, p[m], pc * np.exp(f / tr[m]), "vapour pressure")
    m = (tr < 0.98) & np.isfinite(rho_l)
    c, f = lsq([t[m] ** e for e in DENSITY_EXP], rho_l[m] / rc - 1.0, rc / rho_l[m])
    res["density"] = c
    band(m, rho_l[m], rc * (1 + f), "density")
    m = (tr < 0.95) & np.isfinite(cp_l)
    c, f = lsq([np.ones(m.sum()), 1 / t[m], t[m], t[m] ** 2], cp_l[m])
    res["cp"] = c
    band(m, cp_l[m], f, "specific heat (Tr < 0.95)")
    lat = h_v - h_l
    m = (tr < 0.98) & np.isfinite(lat)
    c, f = lsq([t[m] ** 0.38, t[m] ** 1.38, t[m] ** 2.38, t[m] ** 3.38], lat[m])
    res["latent"] = c
    band(m, lat[m], f, "latent heat")
    m = (tr < 0.98) & np.isfinite(sig) & (sig > 0)
    c, f = lsq([t[m] ** 1.26, t[m] ** 2.26, t[m] ** 3.26, t[m] ** 4.26], sig[m])
    res["sigma"] = c
    band(m, sig[m], f, "surface tension")
    m = (tr < 0.95) & np.isfinite(mu_l)
    c, f = lsq([tr[m] ** e for e in VISCOSITY_EXP], np.log(mu_l[m]), np.ones(m.sum()))
    res["viscosity"] = c
    band(m, mu_l[m], np.exp(f), "viscosity (Tr < 0.95)")
    return res, report


def wagner(r, t_k):
    tr = t_k / r["tc"]
    t = 1.0 - tr
    a, b, c, d = r["psat"]
    return r["pc"] * math.exp((a * t + b * t**1.5 + c * t**2.5 + d * t**5) / tr)


def rust(name, r, rep):
    mm = FLUIDS[name][1]
    omega = -math.log10(wagner(r, 0.7 * r["tc"]) / r["pc"]) - 1.0
    lo, hi = r["t_min"], r["tc"]
    for _ in range(80):  # normal boiling point by bisection
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if wagner(r, mid) < 101325.0 else (lo, mid)
    vc = mm / r["rho_c"]
    zc = r["pc"] * vc / (R * r["tc"])
    f = lambda v: ", ".join(f"{x:.10e}" for x in v)
    const = name.upper().replace("-", "_")
    src = "the NIST table" if FLUIDS[name][0] else "the table of isooctane_table.py"
    out = [f"/// {name}: largest relative errors of the fits against {src}:"]
    out += [f"/// {e}." for e in rep]
    out.append(f"pub const {const}: LiquidData = LiquidData {{")
    out.append(f'    name: "{name}",')
    out.append(f"    molar_mass: {mm!r},")
    out.append(f"    tc: {float(r["tc"])!r},")
    out.append(f"    pc: {float(r["pc"])!r},")
    out.append(f"    rho_c: {r['rho_c']:.6f},")
    out.append(f"    vc: {vc:.6e},")
    out.append(f"    zc: {zc:.6f},")
    out.append(f"    omega: {omega:.6f},")
    out.append(f"    tb: {0.5 * (lo + hi):.4f},")
    out.append(f"    t_min: {float(r["t_min"])!r},")
    for k in ["psat", "density", "cp", "latent", "sigma", "viscosity"]:
        out.append(f"    {k}: [{f(r[k])}],")
    out.append("};")
    return "\n".join(out)


def main(argv):
    if argv[1:] == ["--urls"]:
        for n, (i, _, lo, hi) in FLUIDS.items():
            if i is None:
                print(n, "python3 tools/data/isooctane_table.py DIR")
            else:
                print(n, URL.format(id=i, lo=lo, hi=hi))
        return
    print("// Written by tools/data/liquid_properties.py from the NIST WebBook")
    print("// saturation tables (SRD 69) and, for iso-octane, isooctane_table.py;")
    print("// do not edit by hand.")
    for n in FLUIDS:
        r, rep = fit(n, f"{argv[1]}/{n}.txt")
        print(rust(n, r, rep))


if __name__ == "__main__":
    main(sys.argv)
