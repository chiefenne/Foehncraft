"""The ideal-gas heat capacity of n-dodecane and its NASA-7 polynomials
(Spray A, SA1, D311), from published data:

- Scott (1974) correlated the ideal-gas properties of the alkanes so that
  each CH2 group adds the same amount: his recommended c_p of n-nonane and
  n-decane (NIST Chemistry WebBook, gas phase thermochemistry, "Scott
  D.W., 1974, 2": D. W. Scott, J. Chem. Phys. 60 (1974) 3144, and the
  Bureau of Mines bulletin 666 of the same year) differ by the same CH2
  increment as n-octane and n-nonane at every temperature, so n-dodecane
  is n-decane plus two increments: c_p(C12) = c_p(C10) + 2 (c_p(C10) -
  c_p(C9)).
- Check: the c_p of the reference equation of state of Lemmon and Huber
  (2004) at 1 Pa (NIST WebBook fluid tables, the ideal gas within 1e-5)
  agrees with the table within 0.11 % at 300, 400, 500, 600 and 700 K, and
  with the fit below within 0.24 % from 300 to 700 K.
- NASA-7 fit as tools/data/vapour_species.py (D320): least squares in the
  relative error, two ranges joined at 1000 K with c_p and dc_p/dT
  continuous, from 273.15 K (the 200 K value of the table is left out: it
  doubles the largest error of the fit and no spray state is that cold).
  a6 from the enthalpy of formation of the gas at 298.15 K, -290.9 kJ/mol
  (Prosen and Rossini 1945), a7 from the standard entropy 622.50 J/(mol K)
  (Stull 1969), both NIST WebBook.
- Lennard-Jones parameters from the critical point and the acentric factor
  of the kernel's fuel data (Ambrose and Walton 1989) by Tee, Gotoh and
  Stewart (1966), as D320.

    python3 tools/data/dodecane_ideal_gas.py            # Rust constant
    python3 tools/data/dodecane_ideal_gas.py --csv      # T, c_p table

The Rust constant goes into `fuels::N_DODECANE.ideal_gas` of
crates/foehncraft-kernels/src/realfluid.rs.
"""

import sys

import numpy as np

R = 8.314462618

# (T K, c_p J/(mol K)), NIST WebBook, Scott 1974 (C111842, C124185)
NONANE = [
    (200.0, 161.92),
    (273.15, 196.77),
    (298.15, 210.4),
    (300.0, 211.42),
    (400.0, 268.82),
    (500.0, 321.54),
    (600.0, 366.1),
    (700.0, 403.34),
    (800.0, 433.88),
    (900.0, 459.82),
    (1000.0, 481.58),
    (1100.0, 499.99),
    (1200.0, 516.31),
    (1300.0, 531.37),
    (1400.0, 543.92),
    (1500.0, 556.47),
]
DECANE = [
    (200.0, 179.08),
    (273.15, 217.9),
    (298.15, 233.1),
    (300.0, 234.18),
    (400.0, 297.98),
    (500.0, 356.43),
    (600.0, 405.85),
    (700.0, 446.43),
    (800.0, 479.9),
    (900.0, 508.36),
    (1000.0, 531.79),
    (1100.0, 551.87),
    (1200.0, 569.44),
    (1300.0, 585.76),
    (1400.0, 598.31),
    (1500.0, 610.86),
]
# c_p of n-dodecane at 1 Pa, Lemmon and Huber (2004), NIST WebBook fluid tables
LEMMON_HUBER = [
    (300.0, 279.710),
    (350.0, 318.296),
    (400.0, 356.369),
    (450.0, 392.536),
    (500.0, 426.124),
    (550.0, 456.882),
    (600.0, 484.817),
    (650.0, 510.085),
    (700.0, 532.916),
]
DFH = -290.9e3
S298 = 622.50
T_LOW, T_MID, T_HIGH = 273.15, 1000.0, 1500.0
TC, PC, OMEGA, M = 658.65, 1.83e6, 0.573, 0.170334


def table():
    return [(t, c10 + 2.0 * (c10 - c9)) for (t, c10), (_, c9) in zip(DECANE, NONANE)]


def fit(t, cp_r):
    """Two c_p/R quartics, continuous with their slope at T_MID."""
    lo, hi = t <= T_MID, t >= T_MID

    def rows(tt):
        return np.column_stack([tt**k for k in range(5)])

    a = np.zeros((lo.sum() + hi.sum(), 10))
    a[: lo.sum(), :5] = rows(t[lo]) / cp_r[lo][:, None]
    a[lo.sum() :, 5:] = rows(t[hi]) / cp_r[hi][:, None]
    b = np.ones(a.shape[0])
    c = np.zeros((2, 10))
    c[0, :5] = [T_MID**k for k in range(5)]
    c[0, 5:] = -c[0, :5]
    c[1, :5] = [k * T_MID ** (k - 1) if k else 0.0 for k in range(5)]
    c[1, 5:] = -c[1, :5]
    kkt = np.block([[2 * a.T @ a, c.T], [c, np.zeros((2, 2))]])
    x = np.linalg.solve(kkt, np.concatenate([2 * a.T @ b, np.zeros(2)]))[:10]
    return x[:5], x[5:], np.max(np.abs(a @ x - 1.0))


def h_rt(a, t):
    return a[0] + a[1] * t / 2 + a[2] * t**2 / 3 + a[3] * t**3 / 4 + a[4] * t**4 / 5


def s_r(a, t):
    return (
        a[0] * np.log(t)
        + a[1] * t
        + a[2] * t**2 / 2
        + a[3] * t**3 / 3
        + a[4] * t**4 / 4
    )


def cp(a, t):
    return R * sum(a[k] * t**k for k in range(5))


def main():
    tab = [(t, c) for t, c in table() if t >= T_LOW]
    if "--csv" in sys.argv:
        for t, c in tab:
            print(f"{t},{c:.2f}")
        return
    t = np.array([x for x, _ in tab])
    c = np.array([y for _, y in tab])
    lo, hi, err = fit(t, c / R)
    lo6 = DFH / R - 298.15 * h_rt(lo, 298.15)
    hi6 = T_MID * (h_rt(lo, T_MID) - h_rt(hi, T_MID)) + lo6
    lo7 = S298 / R - s_r(lo, 298.15)
    hi7 = s_r(lo, T_MID) + lo7 - s_r(hi, T_MID)
    dev = max(abs(cp(lo, tt) / v - 1.0) for tt, v in LEMMON_HUBER)
    sigma = (2.3551 - 0.087 * OMEGA) * (TC / (PC / 101325.0)) ** (1 / 3)
    eps = (0.7915 + 0.1693 * OMEGA) * TC
    f = lambda v: ", ".join(f"{x:.10e}" for x in v)
    print(
        f"// largest c_p error of the fit {err:.1e} against the table, "
        f"{dev:.1e} against Lemmon and Huber (2004), 300 to 700 K"
    )
    print("ideal_gas: Species {")
    print(f"    molar_mass: {M},")
    print(f"    t_range: [{T_LOW}, {T_MID}, {T_HIGH}],")
    print(f"    low: [{f(list(lo) + [lo6, lo7])}],")
    print(f"    high: [{f(list(hi) + [hi6, hi7])}],")
    print(f"    well_depth: {eps:.2f},")
    print(f"    diameter: {sigma:.4f},")
    print("    dipole: 0.0,")
    print("    geometry: 2,")
    print("},")


if __name__ == "__main__":
    main()
