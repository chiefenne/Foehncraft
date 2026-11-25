"""Write the species data of the vapours of the spray liquids (D320),
crates/foehncraft-flow/data/species-vapours.toml, from published data:

- NASA-7 polynomials fitted here to the ideal-gas heat capacities of the
  NIST Chemistry WebBook (gas phase thermochemistry, Scott 1974 for the
  alkanes), two ranges joined at 1000 K with c_p and dc_p/dT continuous;
  a6 from the enthalpy of formation of the gas at 298.15 K (NIST WebBook).
  The entropy constant a7 is set from the standard entropy when NIST lists
  one, otherwise 0: such a species must not be used for chemical
  equilibrium (the record says so in `note`).
- Lennard-Jones parameters from the critical constants and the acentric
  factor with the correlation of Tee, Gotoh and Stewart (1966):
  sigma = (2.3551 - 0.087 omega) (Tc/pc)^(1/3) (Angstrom, pc in atm),
  epsilon/k = (0.7915 + 0.1693 omega) Tc.

    python3 tools/data/vapour_species.py DIR > crates/foehncraft-flow/data/species-vapours.toml

DIR holds <file>_cpgas.csv (T K, c_p J/(mol K)) of each species.
"""

import sys

import numpy as np

R = 8.314462618
WEIGHTS = {"H": 1.00794, "C": 12.0107, "O": 15.9994, "N": 14.0067}

# name: (file, composition, dfH298 J/mol, S298 J/(mol K) or None, Tc K, pc Pa, omega, source)
SPECIES = {
    "NC7H16": ("heptane", {"C": 7, "H": 16}, -187.8e3, None, 540.2, 2.73576e6, 0.349216,  # omega of the liquid fit (D320)
               "n-heptane: Cp,gas Scott 1974, dfH Prosen and Rossini 1945 (NIST WebBook); Tc, pc NIST fluid tables"),
    "IC8H18": ("isooctane", {"C": 8, "H": 18}, -224.1e3, None, 543.9, 2.57e6,
               ("antoine", 3.93679, 1257.84, -52.415),
               "iso-octane: Cp,gas Scott 1974, dfH Prosen and Rossini 1945, Tc and pc NIST WebBook; omega from the NIST Antoine data at 0.7 Tc"),
}
T_MID = 1000.0


def fit(t, cp_r):
    """Two NASA-7 c_p/R polynomials, continuous with their slope at T_MID."""
    lo, hi = t <= T_MID, t >= T_MID
    def rows(tt):
        return np.column_stack([tt**k for k in range(5)])
    a = np.zeros((len(t[lo]) + len(t[hi]), 10))
    a[: lo.sum(), :5] = rows(t[lo]) / cp_r[lo][:, None]
    a[lo.sum():, 5:] = rows(t[hi]) / cp_r[hi][:, None]
    b = np.ones(a.shape[0])
    c = np.zeros((2, 10))
    c[0, :5] = [T_MID**k for k in range(5)]
    c[0, 5:] = -c[0, :5]
    c[1, :5] = [k * T_MID ** (k - 1) if k else 0.0 for k in range(5)]
    c[1, 5:] = -c[1, :5]
    kkt = np.block([[2 * a.T @ a, c.T], [c, np.zeros((2, 2))]])
    x = np.linalg.solve(kkt, np.concatenate([2 * a.T @ b, np.zeros(2)]))[:10]
    return x[:5], x[5:], np.max(np.abs(a @ x - 1.0))


def h_rt_poly(a, t):
    return a[0] + a[1] * t / 2 + a[2] * t**2 / 3 + a[3] * t**3 / 4 + a[4] * t**4 / 5


def s_r_poly(a, t):
    return a[0] * np.log(t) + a[1] * t + a[2] * t**2 / 2 + a[3] * t**3 / 3 + a[4] * t**4 / 4


def acentric(om, tc, pc):
    """The acentric factor, or from Antoine data log10(p/bar) = A - B/(T + C) at 0.7 Tc."""
    if not isinstance(om, tuple):
        return om
    _, a, b, c = om
    p = 1e5 * 10 ** (a - b / (0.7 * tc + c))
    return -np.log10(p / pc) - 1.0


def single(t, cp_r):
    """One c_p/R polynomial over the whole range (both NASA-7 ranges equal)."""
    a = np.column_stack([t**k for k in range(5)]) / cp_r[:, None]
    x, *_ = np.linalg.lstsq(a, np.ones(len(t)), rcond=None)
    return x, x.copy(), np.max(np.abs(a @ x - 1.0))


def main(d):
    print("# Vapours of the spray liquids (D320): NASA-7 fitted to the NIST WebBook")
    print("# ideal-gas heat capacities, Lennard-Jones data of Tee, Gotoh and Stewart")
    print("# (1966). Written by tools/data/vapour_species.py; do not edit by hand.")
    print("# Units as in species-gri30.toml.")
    for name, (f, comp, dfh, s298, tc, pc, om, src) in SPECIES.items():
        data = np.loadtxt(f"{d}/{f}_cpgas.csv", delimiter=",")
        t, cp_r = data[:, 0], data[:, 1] / R
        lo, hi, err = single(t, cp_r)
        if err > 2e-3:
            lo, hi, err = fit(t, cp_r)
        om = acentric(om, tc, pc)
        lo6 = dfh / R - 298.15 * h_rt_poly(lo, 298.15)
        hi6 = T_MID * (h_rt_poly(lo, T_MID) - h_rt_poly(hi, T_MID)) + lo6
        lo7 = (s298 / R - s_r_poly(lo, 298.15)) if s298 is not None else 0.0
        hi7 = s_r_poly(lo, T_MID) + lo7 - s_r_poly(hi, T_MID) if s298 is not None else 0.0
        m = sum(WEIGHTS[e] * n for e, n in comp.items()) / 1000.0
        sigma = (2.3551 - 0.087 * om) * (tc / (pc / 101325.0)) ** (1 / 3)
        eps = (0.7915 + 0.1693 * om) * tc
        f7 = lambda v: ", ".join(f"{x:.10e}" for x in v)
        print()
        print(f"# {src}; largest c_p error of the fit {err:.1e}; acentric factor {om:.4f}")
        if s298 is None:
            print("# entropy constant a7 not set (no standard entropy in the source): not for chemical equilibrium")
        print(f'[species."{name}"]')
        print("composition = { " + ", ".join(f"{e} = {float(n)!r}" for e, n in comp.items()) + " }")
        print(f"molar_mass = {m:.8f}")
        print(f"temperature = [200.0, {T_MID}, 1500.0]")
        print("coefficients = [")
        print(f"  [{f7(list(lo) + [lo6, lo7])}],")
        print(f"  [{f7(list(hi) + [hi6, hi7])}],")
        print("]")
        print('geometry = "nonlinear"')
        print(f"well_depth = {eps:.2f}")
        print(f"diameter = {sigma:.4f}")



if __name__ == "__main__":
    main(sys.argv[1])
