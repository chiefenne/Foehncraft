"""The conservation ledger as data (gallery entry of PR 238, D119): the
terms of the kinetic-energy and the thermal-energy account of a run, read
from the `balance.last.ledger` of its result.json, and the closure of every
account (imbalance over the largest term).

    pip install matplotlib
    python3 -m tools.plot.ledger_chart OUT_DIR RUN_DIR

e.g.

    target/release/foehn flow examples/cases/flow/cavity_heated.toml --out output/cavity_heated
    python3 -m tools.plot.ledger_chart docs/figures/gallery output/cavity_heated

The figure is ledger_cavity.png. Zero terms are left out. The imbalance of
the kinetic-energy account is the numerical dissipation of the scheme (D90),
that of mass and thermal energy is round-off. Sources are positive, sinks negative; the
storage is drawn with the sign it has on the right-hand side of the
closure (storage = sum of the terms).
"""
import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

KIND_COLOUR = {
    "boundary": "#2f6fd0",
    "source": "#e08a1e",
    "dissipation": "#c8372d",
    "exchange": "#2e9a5a",
    "storage": "#7a7a7a",
}


def ledger(run_dir):
    """The accounts of the last balance of the last solve of a run."""
    with open(os.path.join(run_dir, "result.json")) as f:
        result = json.load(f)
    solve = result["solves"][-1]
    return solve["balance"]["last"]["ledger"]


def rows(account):
    """(label, value, kind) of the closed terms of an account and its storage."""
    out = []
    for t in account["terms"]:
        if t["kind"] == "part":
            continue
        if t["value"] == 0.0:
            continue
        label = t["name"] if t.get("boundary") is None else f"{t['name']} ({t['boundary']})"
        out.append((label, t["value"], t["kind"]))
    out.append(("storage", account["storage"], "storage"))
    return out


def terms_panel(ax, account, title):
    data = rows(account)
    labels = [d[0] for d in data][::-1]
    values = [d[1] for d in data][::-1]
    colours = [KIND_COLOUR[d[2]] for d in data][::-1]
    ax.barh(range(len(values)), values, color=colours)
    ax.set_yticks(range(len(values)))
    ax.set_yticklabels(labels, fontsize=7)
    ax.axvline(0, color="#444", lw=0.6)
    ax.set_title(title, fontsize=9)
    ax.set_xlabel(f"{account['unit']} per second", fontsize=8)
    ax.tick_params(axis="x", labelsize=7)


def closure_panel(ax, accounts):
    names = [f"{a['quantity']}" for a in accounts][::-1]
    floor = 1e-16
    vals = [max(abs(a["relative"]), floor) for a in accounts][::-1]
    ax.barh(range(len(vals)), vals, color="#8a4fc8")
    ax.set_xscale("log")
    ax.set_xlim(floor, 1.0)
    ax.set_yticks(range(len(vals)))
    ax.set_yticklabels(names, fontsize=8)
    ax.set_title("closure: |imbalance| over the largest term", fontsize=9)
    ax.set_xlabel("kinetic energy: the numerical dissipation (D90)", fontsize=7)
    ax.tick_params(axis="x", labelsize=7)
    ax.axvline(1e-6, color="#888", lw=0.6, ls="--")


def main(out_dir, run_dir):
    accounts = ledger(run_dir)
    by = {a["quantity"]: a for a in accounts}
    fig, axs = plt.subplots(1, 3, figsize=(12.1, 3.6), gridspec_kw={"width_ratios": [1, 1, 0.8]})
    terms_panel(axs[0], by["kinetic energy"], "kinetic energy: the lid drives, viscosity dissipates")
    terms_panel(axs[1], by["thermal energy"], "thermal energy: the heat balance as rows")
    closure_panel(axs[2], accounts)
    fig.tight_layout()
    os.makedirs(out_dir, exist_ok=True)
    fig.savefig(os.path.join(out_dir, "ledger_cavity.png"), dpi=100)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
