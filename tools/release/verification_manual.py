#!/usr/bin/env python3
"""The verification manual (D68), generated from the verification cases
and their pass bands.

    python3 tools/release/verification_manual.py [OUTPUT] [--docs DOCS]

OUTPUT is the output directory of a catalogue run (default `output`, as
`python3 -m tools.check.ci_cases` writes it: OUTPUT/<case>/ per case file). The script
reads every case file of examples/cases/{flow,scalar,poisson} (the GPU
variants aside), its pass bands (`[verification]`) and the report of its
run, and writes

- DOCS/verification/manual.md: one table per family of cases (the pass
  band of every check and the value of the run), the discretisation
  uncertainty of the sweeps (D67) and the cases without a run;
- DOCS/verification/results.json: the same as data (format
  `foehncraft.verification-manual`, version 1.0);
- DOCS/figures/verification/orders.svg: the observed orders against their
  bands;
- with the default DOCS, docs/verification/status.md and grades.json: the
  grades of the features (D199, tools/release/verification_status.py).

DOCS defaults to `docs`. Regenerate it for every release and after a
change of a scheme, from a run of all catalogue sets:

    for s in linux-1 linux-2 linux-3 long; do python3 -m tools.check.ci_cases $s target/release/foehn; done
    python3 tools/release/verification_manual.py output

The output depends only on the reports (no dates, no timings), so the
same run gives the same files.
"""
from __future__ import annotations

import json
import pathlib
import re
import sys
import tomllib

ROOT = pathlib.Path(__file__).resolve().parents[2]
FAMILIES = [
    ("flow", "result.json", "Flow cases (`foehn flow`)"),
    ("scalar", "scalar.json", "Scalar transport (`foehn scalar`)"),
    ("poisson", "poisson.json", "Pressure equation (`foehn poisson`)"),
]
FORMAT, VERSION = "foehncraft.verification-manual", "1.0"


def describe(path: pathlib.Path) -> str:
    """The first sentence of the leading comment of a case file, without
    its stage note."""
    words = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.startswith("#"):
            break
        text = line.lstrip("#").strip()
        if not text:
            if words:
                break
            continue
        words.append(text)
    text = " ".join(words).split(". ")[0]
    head = text.split(": ")[0]
    if len(head) >= 30:
        text = head
    text = re.sub(r"[;,]\s*stage [^)]*\)", ")", text)
    text = re.sub(r"\s*\(stage [^)]*\)", "", text)
    return text.rstrip(".:") or path.stem


def band(lo, hi) -> str:
    return f"{lo:g} – {hi:g}"


def check(name, band_text, value, ok, unit="") -> dict:
    return {"check": name, "band": band_text, "value": value, "unit": unit, "passed": bool(ok)}


def flow_checks(case: dict, rep: dict) -> list[dict]:
    """The checks of a flow case, as the solver makes them (case.rs,
    unsteady.rs), with the band and the value of the run."""
    v = case.get("verification", {})
    out = []
    solves = rep.get("solves", [])
    if v.get("expected_order") and rep.get("orders"):
        lo, hi = v["expected_order"]
        o = rep["orders"][-1]
        out.append(check("observed order in space (finest pair)", band(lo, hi), o, lo <= o <= hi))
    if v.get("expected_time_order") and rep.get("time_orders"):
        lo, hi = v["expected_time_order"]
        o = rep["time_orders"][-1]
        out.append(check("observed order in time (finest pair)", band(lo, hi), o, lo <= o <= hi))
    if v.get("scalar_order"):
        lo, hi = v["scalar_order"]
        for name, orders in rep.get("scalar_orders", []):
            if orders:
                o = orders[-1]
                what = "in time" if name.endswith(" (time)") else "in space"
                out.append(check(f"observed order of scalar {name.removesuffix(' (time)')} {what} (finest pair)",
                                 band(lo, hi), o, lo <= o <= hi))
    if v.get("scalar_max_error") is not None:
        for name in {sc["name"] for s in solves[-1:] for sc in s.get("scalars", [])}:
            l2 = [sc["l2"] for sc in solves[-1].get("scalars", []) if sc["name"] == name and sc.get("l2") is not None]
            if l2:
                out.append(check(f"L2 error of scalar {name} (finest mesh)", f"≤ {v['scalar_max_error']:g}",
                                 l2[0], l2[0] <= v["scalar_max_error"]))
    if v.get("max_error") is not None:
        errs = [s["errors"]["l2"] for s in solves if s.get("errors")]
        if errs:
            e = max(errs)
            out.append(check("velocity error L2 (largest of the levels)", f"≤ {v['max_error']:g}",
                             e, e <= v["max_error"], "m/s"))
    if v.get("reference_tolerance") is not None:
        devs = [max(s["reference"]["u_max_dev"], s["reference"]["v_max_dev"])
                for s in solves if s.get("reference")]
        if devs:
            ref = next(s["reference"]["reference"] for s in solves if s.get("reference"))
            d = max(devs)
            out.append(check(f"largest deviation from {ref}", f"≤ {v['reference_tolerance']:g}",
                             d, d <= v["reference_tolerance"], "of U_lid"))
    for s in solves[-1:]:
        b = s.get("benchmark")
        if not b:
            continue
        tol = v.get("benchmark_tolerance", 0.0)
        for x in b.get("values", []):
            r = x["range"]
            out.append(check(f"{r['name']} ({b['benchmark']})", band(r["lower"], r["upper"]),
                             x["value"], x["deviation"] <= tol))
    if v.get("require_convergence") and solves:
        steady = [s for s in solves if not s.get("transient")]
        if steady:
            ok = all(s.get("converged") for s in steady)
            out.append(check("every steady solve converged", "converged",
                             "yes" if ok else "no", ok))
    if solves:
        m = max(abs(s.get("mass_imbalance", 0.0)) for s in solves)
        out.append(check("net mass flow out of the domain", "≤ 1e-06", m, m <= 1e-6,
                         "relative"))
    return out


def scalar_checks(case: dict, rep: dict) -> list[dict]:
    v = case.get("verification", {})
    out = []
    if v.get("expected_order") and rep.get("orders"):
        lo, hi = v["expected_order"]
        o = rep["orders"][-1]
        out.append(check("observed order (finest pair)", band(lo, hi), o, lo <= o <= hi))
    bounded = [("", v["bounds"])] if v.get("bounds") else []
    bounded += [(x.get("name", ""), x["bounds"]) for x in case.get("variants", []) if x.get("bounds")]
    for name, (lo, hi) in bounded:
        r = [s["range"] for s in rep.get("solves", [])
             if s.get("range") and (not name or s.get("variant") == name)]
        if r:
            a, b = min(x[0] for x in r), max(x[1] for x in r)
            what = f"range of the solution ({name})" if name else "range of the solution"
            out.append(check(what, f"{band(lo, hi)} (± 1e-6)", f"{fmt(a)} … {fmt(b)}",
                             lo - 1e-6 <= a and b <= hi + 1e-6))
    solves = rep.get("solves", [])
    if solves:
        ok = all(s.get("outcome", {}).get("converged") for s in solves)
        out.append(check("every solve converged", "converged", "yes" if ok else "no", ok))
    return out


def poisson_checks(case: dict, rep: dict) -> list[dict]:
    mx = case.get("study", {}).get("max_iteration_ratio")
    out = []
    for solver, ratio in rep.get("iteration_ratios", []):
        out.append(check(f"{solver}: iterations finest / coarsest level",
                         f"≤ {mx:g}" if mx else "reported", ratio, mx is None or ratio <= mx))
    solves = rep.get("solves", [])
    if solves:
        d = max((s.get("diff_to_reference") or 0.0) for s in solves)
        tol = case.get("study", {}).get("max_difference", 1e-6)
        out.append(check("largest difference to the reference solution", f"≤ {tol:g}", d, d <= tol,
                         "relative"))
    return out


CHECKS = {"flow": flow_checks, "scalar": scalar_checks, "poisson": poisson_checks}


def meshes(family: str, rep: dict) -> str:
    solves = rep.get("solves", [])
    if family == "flow":
        cells = sorted({s.get("cells", 0) for s in solves})
    else:
        cells = sorted({s.get("cells", 0) for s in solves if s.get("cells")})
    if not cells:
        return "-"
    return " … ".join(f"{c:,}".replace(",", " ") for c in (cells[0], cells[-1])) if len(cells) > 1 \
        else f"{cells[0]:,}".replace(",", " ")


def collect(output: pathlib.Path) -> dict:
    families = []
    for family, report, title in FAMILIES:
        cases = []
        for path in sorted((ROOT / "examples" / "cases" / family).glob("*.toml")):
            if path.stem.endswith("_gpu"):
                continue
            case = tomllib.loads(path.read_text(encoding="utf-8"))
            rp = output / path.stem / report
            entry = {
                "case": path.stem,
                "file": str(path.relative_to(ROOT)),
                "description": describe(path),
            }
            if not rp.is_file():
                entry["status"] = "not run"
                cases.append(entry)
                continue
            rep = json.loads(rp.read_text(encoding="utf-8"))
            failures = rep.get("failures") or [
                d.get("message", "") for d in rep.get("diagnostics", []) if d.get("severity") == "error"
            ]
            status = rep.get("status") or ("failed" if failures else "passed")
            checks = CHECKS[family](case, rep)
            passed = status in ("finished", "passed") and not failures and all(c["passed"] for c in checks)
            entry.update(
                status="passed" if passed else "failed",
                cells=meshes(family, rep),
                checks=checks,
                orders=rep.get("orders", []),
                time_orders=rep.get("time_orders", []),
                uncertainty=rep.get("uncertainty", []),
                failures=failures,
            )
            if family == "flow":
                prov = rep.get("provenance", {})
                entry["provenance"] = {k: prov.get(k) for k in ("version", "commit", "materials_database")}
            cases.append(entry)
        families.append({"family": family, "title": title, "cases": cases})
    versions = sorted({
        (c["provenance"].get("version"), c["provenance"].get("commit"))
        for f in families for c in f["cases"] if c.get("provenance")
    }, key=str)
    return {"format": FORMAT, "format_version": VERSION, "runs_of": [
        {"version": a, "commit": b} for a, b in versions], "families": families}


def fmt(x) -> str:
    if isinstance(x, str):
        return x
    if isinstance(x, float):
        if x == 0.0:
            return "0"
        if 1e-3 <= abs(x) < 1e4:
            return f"{x:.4g}"
        return f"{x:.3e}"
    return str(x)


def markdown(data: dict) -> str:
    out = []
    total = [c for f in data["families"] for c in f["cases"]]
    ran = [c for c in total if c["status"] != "not run"]
    passed = [c for c in ran if c["status"] == "passed"]
    runs = ", ".join(
        f"{r['version']} ({(r['commit'] or 'no commit')[:12]})" for r in data["runs_of"]
    ) or "unknown"
    out += [
        "# Verification manual",
        "",
        "<!-- generated by tools/release/verification_manual.py (D68); do not edit -->",
        "",
        "Every verification case of the [case catalogue](../cases.md) with its pass "
        "bands and the result of its run. The cases check themselves (a run outside "
        "its band exits with code 1, R5); this page collects the bands and the values. "
        "The explanations stay on the stage pages ([stage 2](stage2.md), "
        "[stage 3](stage3.md), [stage 4](stage4.md), [stage 5](stage5.md), [stage 6](stage6.md)). "
        "The same data as JSON: [results.json](results.json).",
        "",
        f"Runs of Foehncraft {runs}: {len(passed)} of {len(ran)} cases passed, "
        f"{len(total) - len(ran)} without a run.",
        "",
        "![Observed orders against their bands](../figures/verification/orders.svg)",
        "",
    ]
    for f in data["families"]:
        out += [f"## {f['title']}", ""]
        out += ["| Case | What it checks | Cells | Check | Band | Value | |", "|---|---|---|---|---|---|---|"]
        for c in f["cases"]:
            if c["status"] == "not run":
                continue
            rows = c["checks"] or [check("the run itself", "exit code 0", c["status"], c["status"] == "passed")]
            for i, k in enumerate(rows):
                head = (f"`{c['case']}`", c["description"], c["cells"]) if i == 0 else ("", "", "")
                value = fmt(k["value"]) + (f" {k['unit']}" if k["unit"] and not isinstance(k["value"], str) else "")
                out.append(f"| {head[0]} | {head[1]} | {head[2]} | {k['check']} | {k['band']} | {value} | "
                           f"{'passed' if k['passed'] else '**failed**'} |")
        out.append("")
        failed = [c for c in f["cases"] if c["status"] == "failed"]
        for c in failed:
            out.append(f"- `{c['case']}` failed: " + ("; ".join(c["failures"]) or "a check above"))
        if failed:
            out.append("")
        sweeps = [c for c in f["cases"] if c.get("uncertainty")]
        for c in sweeps:
            out += [
                f"### Discretisation uncertainty of `{c['case']}`",
                "",
                "Three finest solutions, Celik et al. (2008) (D67): observed order p, "
                "the extrapolated value, the band φ1 ± U with U = 1.25 |φ1 − φ2| / (r^p − 1) "
                "and the GCI U / |φ1| where |φ1| is larger than U.",
                "",
                "| Quantity | φ1 (finest) | p | Extrapolated | ± U | GCI | Verdict |",
                "|---|---|---|---|---|---|---|",
            ]
            for e in c["uncertainty"]:
                gci = e.get("gci")
                b = e.get("band")
                show_gci = gci is not None and b is not None and b < abs(e["values"][0])
                out.append(
                    f"| {e['quantity']} | {e['values'][0]:.6e} | "
                    f"{'-' if e.get('order') is None else format(e['order'], '.3f')} | "
                    f"{'-' if e.get('extrapolated') is None else format(e['extrapolated'], '.6e')} | "
                    f"{'-' if b is None else format(b, '.2e')} | "
                    f"{format(100 * gci, '.3f') + ' %' if show_gci else '-'} | {e['verdict']} |"
                )
            out.append("")
    missing = [c for c in total if c["status"] == "not run"]
    if missing:
        out += ["## Cases without a run", "",
                "These case files had no report in the output directory "
                "(for example the `long` set of `tools/check/ci_cases.py`):", ""]
        out += [f"- `{c['file']}`: {c['description']}" for c in missing]
        out.append("")
    return "\n".join(out).rstrip("\n") + "\n"


def svg(data: dict) -> str:
    """Observed orders (finest pair) of the cases with an order band."""
    rows = []
    for f in data["families"]:
        for c in f["cases"]:
            for k in c.get("checks", []):
                if k["check"].startswith("observed order"):
                    lo, hi = (float(x) for x in k["band"].split(" – "))
                    tag = " (time)" if "time" in k["check"] else ""
                    rows.append((c["case"] + tag, lo, hi, k["value"], k["passed"]))
    w, left, right, top, row = 720, 230, 30, 46, 18
    h = top + row * max(len(rows), 1) + 50
    x0, x1 = 1.4, 2.6
    sx = lambda v: left + (min(max(v, x0), x1) - x0) / (x1 - x0) * (w - left - right)  # noqa: E731
    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}" '
        'font-family="Helvetica, Arial, sans-serif" font-size="11">',
        f'<rect width="{w}" height="{h}" fill="white"/>',
        f'<text x="{w / 2}" y="20" text-anchor="middle" font-size="13" font-weight="bold">'
        "Observed orders of the verification cases against their bands</text>",
    ]
    for v in (1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6):
        out.append(f'<line x1="{sx(v):.1f}" y1="{top - 6}" x2="{sx(v):.1f}" y2="{h - 40}" stroke="#ddd"/>')
        out.append(f'<text x="{sx(v):.1f}" y="{h - 26}" text-anchor="middle">{v:g}</text>')
    out.append(f'<text x="{(left + w - right) / 2}" y="{h - 8}" text-anchor="middle">observed order</text>')
    for i, (name, lo, hi, val, ok) in enumerate(rows):
        y = top + row * i + row / 2
        out.append(f'<text x="{left - 8}" y="{y + 4:.1f}" text-anchor="end">{name}</text>')
        out.append(f'<rect x="{sx(lo):.1f}" y="{y - 5:.1f}" width="{sx(hi) - sx(lo):.1f}" height="10" '
                   'fill="#1f6fb4" fill-opacity="0.25"/>')
        colour = "#1f6fb4" if ok else "#c0392b"
        out.append(f'<circle cx="{sx(val):.1f}" cy="{y:.1f}" r="4" fill="{colour}"/>')
    out.append("</svg>")
    return "\n".join(out) + "\n"


def main(argv: list[str]) -> int:
    args = [a for a in argv if not a.startswith("--")]
    docs = ROOT / "docs"
    if "--docs" in argv:
        docs = pathlib.Path(argv[argv.index("--docs") + 1])
        args = [a for a in args if a != argv[argv.index("--docs") + 1]]
    output = pathlib.Path(args[0]) if args else ROOT / "output"
    data = collect(output)
    (docs / "verification").mkdir(parents=True, exist_ok=True)
    (docs / "figures" / "verification").mkdir(parents=True, exist_ok=True)
    (docs / "verification" / "manual.md").write_text(markdown(data), encoding="utf-8")
    (docs / "verification" / "results.json").write_text(
        json.dumps(data, indent=1, ensure_ascii=False, allow_nan=False, default=str) + "\n", encoding="utf-8")
    (docs / "figures" / "verification" / "orders.svg").write_text(svg(data), encoding="utf-8")
    if docs.resolve() == (ROOT / "docs").resolve():
        # the grades of D199 follow the new results (docs/verification/status.md)
        sys.path.insert(0, str(ROOT))
        from tools.release import verification_status  # noqa: PLC0415

        verification_status.write()
    cases = [c for f in data["families"] for c in f["cases"]]
    ran = [c for c in cases if c["status"] != "not run"]
    bad = [c["case"] for c in ran if c["status"] != "passed"]
    print(f"{len(ran) - len(bad)} of {len(ran)} cases passed, {len(cases) - len(ran)} without a run; "
          f"wrote {docs / 'verification' / 'manual.md'}")
    if bad:
        print("failed: " + ", ".join(bad))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
