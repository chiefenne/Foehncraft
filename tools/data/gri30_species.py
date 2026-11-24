"""Write the species data file of SF1 (crates/foehncraft-flow/data/species-gri30.toml)
from the GRI-Mech 3.0 thermo and transport data.

GRI-Mech 3.0 (Smith, Golden, Frenklach, Moriarty, Eiteneer, Goldenberg,
Bowman, Hanson, Song, Gardiner, Lissianski and Qin, 1999,
http://combustion.berkeley.edu/gri-mech/) is free to use with its citation.
The input is the mechanism in the YAML form of D103 (`gri30.yaml`, the
GRI-Mech THERMO and TRANSPORT records converted one to one); only the
species records are read, no reaction data and no code.

    python3 tools/data/gri30_species.py path/to/gri30.yaml > crates/foehncraft-flow/data/species-gri30.toml
"""

import sys

import yaml

# IUPAC standard atomic weights as in GRI-Mech 3.0 (g/mol)
WEIGHTS = {"H": 1.00794, "C": 12.0107, "O": 15.9994, "N": 14.0067, "Ar": 39.948}


def main(path):
    d = yaml.safe_load(open(path))
    out = [
        "# Species of GRI-Mech 3.0 (Smith et al. 1999, http://combustion.berkeley.edu/gri-mech/):",
        "# NASA 7-coefficient polynomials of c_p/R, h/(RT), s/R (two ranges) and the",
        "# Lennard-Jones transport data of the GRI-Mech TRANSPORT file.",
        "# Written by tools/data/gri30_species.py; do not edit by hand.",
        "# Units: molar_mass kg/mol, temperatures K, well_depth (epsilon/k_B) K,",
        "# diameter Angstrom, dipole Debye, polarizability Angstrom^3.",
        "",
    ]
    for s in d["species"]:
        th = s["thermo"]
        assert th["model"] == "NASA7"
        tr = s["transport"]
        m = sum(WEIGHTS[e] * n for e, n in s["composition"].items()) / 1000.0
        out.append(f'[species.{s["name"]}]')
        comp = ", ".join(f"{e} = {n:g}" for e, n in s["composition"].items())
        out.append(f"composition = {{ {comp} }}")
        out.append(f"molar_mass = {round(m, 12)!r}")
        out.append(f"temperature = {th['temperature-ranges']!r}")
        rows = ",\n  ".join("[" + ", ".join(repr(float(x)) for x in r) + "]" for r in th["data"])
        out.append(f"coefficients = [\n  {rows},\n]")
        out.append(f'geometry = "{tr["geometry"]}"')
        out.append(f"well_depth = {float(tr['well-depth'])!r}")
        out.append(f"diameter = {float(tr['diameter'])!r}")
        for k in ("dipole", "polarizability", "rotational-relaxation"):
            if k in tr:
                out.append(f"{k.replace('-', '_')} = {float(tr[k])!r}")
        out.append("")
    sys.stdout.write("\n".join(out))


if __name__ == "__main__":
    main(sys.argv[1])
