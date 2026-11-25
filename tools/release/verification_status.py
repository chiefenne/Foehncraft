#!/usr/bin/env python3
"""The verification status page (D199): a grade for every feature.

    python3 -m tools.release.verification_status           write grades.json
    python3 -m tools.release.verification_status --all     also the page and the chart
    python3 -m tools.release.verification_status --check   fail when grades.json is out of date

The feature table `docs/verification/features.toml` is written by hand and
reviewed in pull requests: one row per switch a user can turn on, the
evidence that checks it (catalogue cases, tests, recorded runs), the known
gaps and the validated range. The grade of a row is computed here from the
table and the last verification manual (`docs/verification/results.json`);
it is never written by hand. The script writes

- docs/verification/status.md: the grades, how each feature was checked
  (with the bands of the manual) and the gaps;
- docs/verification/grades.json: the same as data (format
  `foehncraft.verification-grades`, version 1.0), for the limitations of
  a run (D199 LIM001, LIM002).

`--check` (CI, also on docs-only changes) refuses a table that names a case
file, test or page that does not exist, leaves a catalogue case without a
feature, or disagrees with grades.json.

Pull requests commit the table and grades.json only (grades.json is
compiled into the solver for LIM001, so its tests need it). The page and the
chart carry totals, so two pull requests that each add a row would always
conflict on them: they are regenerated on main after every merge by the
workflow `docs-index.yml` (`--all`), like the number tables of the docs.

Grade rules (D199, also printed on the page):

- an evidence item *holds* when its test exists (tests run on every pull
  request), its recorded run is documented, or its case passed in the last
  manual run or is in a CI case set of `tools/check/ci_cases.py`;
- **unverified**: no item holds, or a case of the feature failed in the
  last manual run;
- **verified**: an `exact` or `manufactured` item holds, a
  `reference-data` item holds where the row is a model (`model = true`),
  no gap is written and no item lacks its run;
- **partial**: everything else (gaps written, only `consistency` checks,
  a model without reference data, a case without a run).

Infrastructure features (restart, logs, run record, Python package, command
line, tool registry, journal and replay, MCP server) have no physics to
compare with, so no exact solution exists for them. A row marked
`infrastructure = true` is **verified** through *contract evidence*
instead: it lists the properties it promises (`contracts = [...]`, from
`PROPERTIES`) and every one needs a holding evidence item of kind
`contract` with that `property` (a test that checks it, run on every pull
request); no gap may be written. Missing contract evidence gives
**partial**. Exact, manufactured and reference-data items still count where
they exist; consistency checks stay what they were.

So a gap written into the table, or a failing case, lowers the grade
automatically; raising it needs evidence.
"""
from __future__ import annotations

import json
import os
import pathlib
import re
import sys
import tomllib

ROOT = pathlib.Path(__file__).resolve().parents[2]
TABLE = ROOT / "docs" / "verification" / "features.toml"
RESULTS = ROOT / "docs" / "verification" / "results.json"
PAGE = ROOT / "docs" / "verification" / "status.md"
GRADES = ROOT / "docs" / "verification" / "grades.json"
FIGURE = ROOT / "docs" / "figures" / "verification" / "grades.svg"
COLOURS = {"verified": "#2e9d5b", "partial": "#e9a23b", "unverified": "#d64541"}
FORMAT, VERSION = "foehncraft.verification-grades", "1.0"
KINDS = ("exact", "manufactured", "reference-data", "consistency", "contract")
# the properties contract evidence can prove (infrastructure rows, D199)
PROPERTIES = {
    "bitwise": "two paths to the same state give bitwise equal results (restart, repeat, backend)",
    "schema": "the output has a fixed, tested format (keys, types, format name and version)",
    "equivalence": "two entry points (command line, Python, tool call, MCP) give the same result",
    "replay": "a recorded session replays to the same files",
    "control": "a control file, signal or port takes effect and leaves a clean state",
    "diagnostics": "errors and refusals are structured diagnostics with a code and a hint",
    "round-trip": "writing and reading back gives the same object",
}
GRADE_ORDER = ("verified", "partial", "unverified")
FAMILIES = ("flow", "scalar", "poisson")
# checks of the manual every flow case has; they say nothing about a feature
GENERIC = ("net mass flow out of the domain", "every steady solve converged", "every solve converged")
RANGE_KEYS = ("y_plus_first", "reynolds", "prandtl", "courant", "skewness", "mach")


class TableError(Exception):
    pass


def ci_sets() -> dict[str, str]:
    """Case name -> the CI set that runs it."""
    sys.path.insert(0, str(ROOT))
    from tools.check import ci_cases  # noqa: PLC0415

    out = {}
    for name in ci_cases.LINUX:
        for _, f in ci_cases.cases(name):
            out.setdefault(pathlib.PurePosixPath(f).stem, name)
    return out


def case_files() -> dict[str, str]:
    """Case name -> repository path, for every catalogue case."""
    out = {}
    for fam in FAMILIES:
        for p in sorted((ROOT / "examples" / "cases" / fam).glob("*.toml")):
            out[p.stem] = p.relative_to(ROOT).as_posix()
    return out


def manual() -> tuple[dict[str, dict], list[dict]]:
    data = json.loads(RESULTS.read_text(encoding="utf-8"))
    cases = {c["case"]: c for f in data["families"] for c in f["cases"]}
    return cases, data.get("runs_of", [])


def test_exists(ref: str) -> bool:
    """`path::function` names a Rust `#[test]` or a pytest `test_` function."""
    path, _, name = ref.partition("::")
    p = ROOT / path
    if not name or not p.is_file():
        return False
    text = p.read_text(encoding="utf-8")
    if p.suffix == ".py":
        return name.startswith("test_") and re.search(rf"^def {re.escape(name)}\(", text, re.MULTILINE) is not None
    pattern = rf"#\[test\](?:\s*#\[[^\]]*\])*\s*(?:pub\s+)?fn {re.escape(name)}\s*\("
    return re.search(pattern, text) is not None


def anchors(path: pathlib.Path) -> set[str]:
    """The anchors of a Markdown page as GitHub makes them: explicit `<a id>`
    and the slugs of the headings (numbered when one repeats)."""
    out, seen = set(), {}
    text = path.read_text(encoding="utf-8")
    out.update(re.findall(r'<a id="([^"]+)"', text))
    fence = False
    for line in text.splitlines():
        if line.startswith("```"):
            fence = not fence
        if fence or not line.startswith("#"):
            continue
        title = line.lstrip("#").strip()
        slug = re.sub(r"[^\w\- ]", "", title.lower()).replace(" ", "-")
        n = seen.get(slug, 0)
        seen[slug] = n + 1
        out.add(slug if n == 0 else f"{slug}-{n}")
    return out


def doc_ok(ref: str) -> bool:
    target, _, anchor = ref.partition("#")
    p = ROOT / target
    return p.is_file() and (not anchor or anchor in anchors(p))


def load_table() -> list[dict]:
    rows = tomllib.loads(TABLE.read_text(encoding="utf-8")).get("feature", [])
    files = case_files()
    seen, errors = set(), []
    for r in rows:
        fid = r.get("id", "?")
        for key in ("id", "area", "title", "evidence"):
            if key not in r:
                errors.append(f"{fid}: missing `{key}`")
        if fid in seen:
            errors.append(f"{fid}: id used twice")
        seen.add(fid)
        if r.get("infrastructure"):
            need = r.get("contracts", [])
            if not need:
                errors.append(f"{fid}: an infrastructure row lists its `contracts`")
            for c in need:
                if c not in PROPERTIES:
                    errors.append(f"{fid}: unknown contract `{c}` (known: {', '.join(PROPERTIES)})")
        elif "contracts" in r:
            errors.append(f"{fid}: `contracts` needs `infrastructure = true`")
        for k in r.get("range", {}):
            if k not in RANGE_KEYS:
                errors.append(f"{fid}: unknown range key `{k}` (known: {', '.join(RANGE_KEYS)})")
        for e in r.get("evidence", []):
            if e.get("kind") not in KINDS:
                errors.append(f"{fid}: evidence kind `{e.get('kind')}` (known: {', '.join(KINDS)})")
            if not ({"case", "test", "run"} & set(e)):
                errors.append(f"{fid}: evidence needs a `case`, `test` or `run`")
            if e.get("kind") == "contract" and e.get("property") not in PROPERTIES:
                errors.append(f"{fid}: contract evidence needs a `property` ({', '.join(PROPERTIES)})")
            elif e.get("kind") != "contract" and "property" in e:
                errors.append(f"{fid}: `property` belongs to contract evidence")
            if "what" not in e:
                errors.append(f"{fid}: evidence without `what`")
            if "case" in e and e["case"] not in files:
                errors.append(f"{fid}: no catalogue case `{e['case']}`")
            if "test" in e and not test_exists(e["test"]):
                errors.append(f"{fid}: no test `{e['test']}` (path::function)")
            if "run" in e and "doc" not in e:
                errors.append(f"{fid}: a recorded run needs the `doc` that records it")
            if "doc" in e and not doc_ok(e["doc"]):
                errors.append(f"{fid}: no page or anchor `{e['doc']}`")
    used = {e["case"] for r in rows for e in r.get("evidence", []) if "case" in e}
    for name in sorted(set(files) - used):
        errors.append(f"catalogue case `{name}` ({files[name]}) checks no feature: add it to a row")
    if errors:
        raise TableError("\n".join(errors))
    return rows


def band_text(entry: dict | None) -> str:
    """The specific checks of a case in the manual, with their bands."""
    if not entry or not entry.get("checks"):
        return ""
    parts = [f"{k['check']} {k['band']}" for k in entry["checks"] if k["check"] not in GENERIC]
    return "; ".join(parts)


def grade(row: dict, results: dict[str, dict], sets: dict[str, str]) -> dict:
    items, auto_gaps, failed = [], [], []
    for e in row["evidence"]:
        state, how = "missing", ""
        case = e.get("case")
        entry = results.get(case) if case else None
        if case and entry and entry.get("status") == "failed":
            state, how = "failed", "failed in the last manual run"
            failed.append(case)
        elif "test" in e:
            state, how = "holds", "test, run on every pull request"
        elif case and entry and entry.get("status") == "passed":
            state, how = "holds", "passed in the last manual run"
        elif case and case in sets:
            state, how = "holds", f"passes in CI (set {sets[case]}), not yet in the manual"
        elif "run" in e:
            state, how = "holds", f"recorded run: {e['run']}"
        if state == "missing" and case:
            auto_gaps.append(f"`{case}` has no run in the verification manual or in CI")
        items.append({**e, "state": state, "how": how, "bands": band_text(entry)})
    holding = {i["kind"] for i in items if i["state"] == "holds"}
    held = {i["property"] for i in items if i["state"] == "holds" and i["kind"] == "contract"}
    lacking = [c for c in row.get("contracts", []) if c not in held] if row.get("infrastructure") else []
    gaps = list(row.get("gaps", [])) + auto_gaps
    if failed or not holding:
        g = "unverified"
    elif row.get("infrastructure"):
        g = "verified" if not lacking and not gaps else "partial"
    elif (holding & {"exact", "manufactured"}) and ("reference-data" in holding or not row.get("model")) \
            and not gaps:
        g = "verified"
    else:
        g = "partial"
    why = []
    if failed:
        why.append("failed in the last manual run: " + ", ".join(f"`{c}`" for c in failed))
    if not holding:
        why.append("no evidence holds")
    if g == "partial":
        if row.get("infrastructure"):
            if lacking:
                why.append("no contract evidence for: " + ", ".join(lacking))
        elif not holding & {"exact", "manufactured"}:
            why.append("no exact or manufactured solution")
        if row.get("model") and "reference-data" not in holding:
            why.append("a model without reference data")
        if gaps:
            why.append("gaps written")
    return {"grade": g, "why": why, "items": items, "gaps": gaps}


def rel(path: str) -> str:
    """A repository path as a link from docs/verification/."""
    target, _, anchor = path.partition("#")
    link = os.path.relpath(ROOT / target, PAGE.parent).replace(os.sep, "/")
    return link + (f"#{anchor}" if anchor else "")


def evidence_md(i: dict, files: dict[str, str]) -> str:
    heads = []
    if "case" in i:
        heads.append(f"[`{i['case']}`]({rel(files[i['case']])})")
    if "test" in i:
        path, _, name = i["test"].partition("::")
        heads.append(f"[`{name}`]({rel(path)})")
    if "run" in i and not heads:
        heads.append(f"recorded run ({i['run']})")
    text = f"{i['kind']}" + (f" ({i['property']})" if i["kind"] == "contract" else "") + ": " + ", ".join(heads) + f": {i['what']}"
    if i["how"].startswith("passes in CI"):
        text += f" ({i['how']})"
    if i["bands"]:
        text += f" (manual: {i['bands']})"
    if "doc" in i:
        text += f" ([results]({rel(i['doc'])}))"
    if i["state"] != "holds":
        text += f" **{i['state']}**"
    return text


def build(rows: list[dict]) -> tuple[str, dict]:
    results, runs = manual()
    sets = ci_sets()
    files = case_files()
    graded = [(r, grade(r, results, sets)) for r in rows]
    count = {g: sum(1 for _, x in graded if x["grade"] == g) for g in GRADE_ORDER}
    run_text = ", ".join(f"{r['version']} ({(r['commit'] or 'no commit')[:12]})" for r in runs) or "unknown"
    out = [
        "# Verification status",
        "",
        "<!-- generated by tools/release/verification_status.py (D199) from features.toml "
        "and results.json; do not edit -->",
        "",
        "Every feature a case file can switch on, with its grade, the evidence that checks it "
        "and the known gaps. The table of features is [features.toml](features.toml), written by "
        "hand and reviewed in pull requests; the grades are computed from it and from the last "
        f"[verification manual](manual.md) (runs of Foehncraft {run_text}), never typed in. "
        "The same as data: [grades.json](grades.json).",
        "",
        f"**{len(graded)} features: {count['verified']} verified, {count['partial']} partial, "
        f"{count['unverified']} unverified.**",
        "",
        "## How a grade is computed",
        "",
        "Kinds of evidence: `exact` (an exact solution, often with its observed order), "
        "`manufactured` (a manufactured solution with its order), `reference-data` "
        "(experiments, DNS, published benchmarks and correlations) and `consistency` "
        "(conservation, equality of two paths, bitwise equality, input checks) and `contract` "
        "(only for infrastructure features: a tested promise such as bitwise restart, a fixed "
        "output schema, replay to identical files, equality of command line and Python).",
        "",
        "An item *holds* when its test exists (tests run on every pull request), its recorded run "
        "is documented, or its case passed in the last manual run or runs in a CI case set.",
        "",
        "| Grade | Rule |",
        "|---|---|",
        "| verified | an `exact` or `manufactured` item holds, a `reference-data` item holds where "
        "the feature is a model, no gap is written and every case has a run |",
        "| partial | some items hold, but gaps are written, or only `consistency` checks hold, or "
        "a model has no reference data, or a case has no run |",
        "| verified (infrastructure) | restart, logs, run record, interfaces: every contract the "
        "row lists (`bitwise`, `schema`, `equivalence`, `replay`, `control`, `diagnostics`, "
        "`round-trip`) has a holding `contract` item, and no gap is written |",
        "| unverified | no item holds, or a case of the feature failed in the last manual run |",
        "",
        "A gap written into the table, or a failing case, lowers the grade at once; raising it "
        "needs evidence. *Range* is the validated range of the feature; a run outside it will "
        "be reported as a limitation (D199).",
        "",
    ]
    areas: list[str] = []
    for r, _ in graded:
        if r["area"] not in areas:
            areas.append(r["area"])
    out += ["## Summary", "", "| Area | Feature | Grade |", "|---|---|---|"]
    for a in areas:
        for r, x in graded:
            if r["area"] == a:
                out.append(f"| {a} | [{r['title']}](#{r['id'].replace('.', '-')}) | {x['grade']} |")
    out.append("")
    for a in areas:
        out += [f"## {a}", ""]
        for r, x in graded:
            if r["area"] != a:
                continue
            out += [f'<a id="{r["id"].replace(".", "-")}"></a>', "", f"### {r['title']}", ""]
            meta = [f"**{x['grade']}**" + (f" ({'; '.join(x['why'])})" if x["why"] else "")]
            meta.append(f"id `{r['id']}`")
            if r.get("switch"):
                meta.append(f"switch `{r['switch']}`")
            if r.get("decisions"):
                meta.append("decisions " + ", ".join(r["decisions"]))
            out += [" · ".join(meta), ""]
            if x["items"]:
                out += [f"- {evidence_md(i, files)}" for i in x["items"]] + [""]
            if x["gaps"]:
                out += ["Gaps:", ""] + [f"- {g}" for g in x["gaps"]] + [""]
            if r.get("notes"):
                out += ["Notes:", ""] + [f"- {n}" for n in r["notes"]] + [""]
            if r.get("infrastructure"):
                out += ["Contracts: " + ", ".join(f"`{c}`" for c in r["contracts"]) + ".", ""]
            if r.get("range"):
                rng = ", ".join(f"{k} {v[0]:g} … {v[1]:g}" for k, v in r["range"].items())
                out += [f"Range: {rng}.", ""]
    data = {
        "format": FORMAT,
        "format_version": VERSION,
        "features": [
            {
                "id": r["id"],
                "title": r["title"],
                "grade": x["grade"],
                "why": x["why"],
                "gaps": x["gaps"],
                "range": r.get("range", {}),
            }
            for r, x in graded
        ],
    }
    return "\n".join(out).rstrip("\n") + "\n", data


def svg(rows: list[dict], data: dict) -> str:
    """The grades per area as stacked bars (gallery 170)."""
    grade_of = {f["id"]: f["grade"] for f in data["features"]}
    areas: list[str] = []
    for r in rows:
        if r["area"] not in areas:
            areas.append(r["area"])
    counts = [(a, [sum(1 for r in rows if r["area"] == a and grade_of[r["id"]] == g) for g in GRADE_ORDER])
              for a in areas]
    w, left, right, top, bar, gap = 800, 230, 40, 70, 20, 8
    h = top + len(areas) * (bar + gap) + 50
    most = max(sum(c) for _, c in counts)
    scale = (w - left - right) / most
    total = [sum(c[i] for _, c in counts) for i in range(3)]
    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}" '
        'font-family="Helvetica, Arial, sans-serif" font-size="12">',
        f'<rect width="{w}" height="{h}" fill="white"/>',
        f'<text x="{w / 2}" y="22" text-anchor="middle" font-size="15" font-weight="bold">'
        f"Verification grades of {len(rows)} features (D199)</text>",
    ]
    x = left
    for g, n in zip(GRADE_ORDER, total, strict=True):
        out.append(f'<rect x="{x}" y="36" width="14" height="14" fill="{COLOURS[g]}"/>')
        out.append(f'<text x="{x + 20}" y="48">{g} ({n})</text>')
        x += 150
    for i, (a, c) in enumerate(counts):
        y = top + i * (bar + gap)
        out.append(f'<text x="{left - 10}" y="{y + bar * 0.7:.1f}" text-anchor="end">{a}</text>')
        x0 = left
        for g, n in zip(GRADE_ORDER, c, strict=True):
            if n:
                out.append(f'<rect x="{x0:.1f}" y="{y}" width="{n * scale:.1f}" height="{bar}" '
                           f'fill="{COLOURS[g]}"/>')
                out.append(f'<text x="{x0 + n * scale / 2:.1f}" y="{y + bar * 0.7:.1f}" '
                           f'text-anchor="middle" fill="white" font-weight="bold">{n}</text>')
            x0 += n * scale
    out.append(f'<text x="{(left + w - right) / 2}" y="{h - 16}" text-anchor="middle" fill="#555">'
               "features per area, grades computed from docs/verification/features.toml</text>")
    out.append("</svg>")
    return "\n".join(out) + "\n"


def render() -> tuple[str, str, str]:
    rows = load_table()
    page, data = build(rows)
    return page, json.dumps(data, indent=1, ensure_ascii=False) + "\n", svg(rows, data)


def write(pages: bool = True) -> None:
    """Regenerate the grades and, with `pages`, the page and the chart
    (also called by verification_manual.py)."""
    page, grades, figure = render()
    GRADES.write_text(grades, encoding="utf-8", newline="\n")
    if pages:
        PAGE.write_text(page, encoding="utf-8", newline="\n")
        FIGURE.parent.mkdir(parents=True, exist_ok=True)
        FIGURE.write_text(figure, encoding="utf-8", newline="\n")


def main(argv: list[str]) -> int:
    try:
        page, grades, figure = render()
    except TableError as e:
        print(f"docs/verification/features.toml:\n{e}", file=sys.stderr)
        return 1
    if "--check" in argv:
        if not GRADES.is_file() or GRADES.read_text(encoding="utf-8").replace("\r\n", "\n") != grades:
            print("out of date: docs/verification/grades.json"
                  "\nregenerate with  python3 -m tools.release.verification_status", file=sys.stderr)
            return 1
        print("docs/verification/grades.json matches docs/verification/features.toml"
              " (the page and the chart are regenerated on main)")
        return 0
    pages = "--all" in argv
    write(pages)
    data = json.loads(grades)
    count = {g: sum(1 for f in data["features"] if f["grade"] == g) for g in GRADE_ORDER}
    print(f"{len(data['features'])} features: " + ", ".join(f"{v} {k}" for k, v in count.items())
          + f"; wrote {(PAGE if pages else GRADES).relative_to(ROOT).as_posix()}"
          + ("" if pages else " (the page and the chart follow on main)"))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
