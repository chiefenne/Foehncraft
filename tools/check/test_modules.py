#!/usr/bin/env python3
"""Every integration test file is compiled: crates with few test programs.

    python3 -m tools.check.test_modules

A crate whose manifest sets `autotests = false` builds its integration
tests as one program, `tests/main.rs`, or a few (`tests/main.rs`,
`tests/main2.rs`, ...), with one `mod` line per file of `tests/` in exactly
one of them (few programs instead of one per file: the test build was the
longest step of CI). Cargo then ignores a file of `tests/` that has no
`mod` line, so its tests would silently never run. Fails (exit code 1) when
a file has no `mod` line, has one in two programs, or a `mod` line names no
file.
"""
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
MOD = re.compile(r"^mod ([a-z0-9_]+);", re.MULTILINE)


def main() -> int:
    bad = []
    for manifest in sorted((ROOT / "crates").glob("*/Cargo.toml")):
        if not re.search(r"^autotests\s*=\s*false", manifest.read_text(encoding="utf-8"), re.MULTILINE):
            continue
        tests = manifest.parent / "tests"
        crate = manifest.parent.relative_to(ROOT).as_posix()
        main_rs = tests / "main.rs"
        if not main_rs.is_file():
            bad.append(f"{crate}: autotests = false but no tests/main.rs")
            continue
        programs = sorted(tests.glob("main*.rs"))
        names = {p.stem for p in programs}
        listed = {}
        for program in programs:
            for m in MOD.findall(program.read_text(encoding="utf-8")):
                listed.setdefault(m, []).append(program.name)
        files = {p.stem for p in tests.glob("*.rs")} - names
        bad += [f"{crate}/tests/{m}.rs: no `mod {m};` in tests/main*.rs (never compiled)" for m in sorted(files - set(listed))]
        bad += [f"{crate}/tests/{m}.rs: `mod {m};` in {' and '.join(listed[m])} (compiled twice)" for m in sorted(listed) if len(listed[m]) > 1]
        bad += [f"{crate}/tests/{', '.join(listed[m])}: `mod {m};` but no tests/{m}.rs" for m in sorted(set(listed) - files)]
    for line in bad:
        print(line, file=sys.stderr)
    if not bad:
        print("every integration test file is a module of one test program of its crate")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
