"""Barcode figure for the README: the longest persistence bars of the intact radical-pair compass and of the polar
control in one simulated bird, under Z2 and Z3 coefficients, against the counting threshold.
Usage: python3 -I scripts/make_figures.py [SEED]   (reads results/confirmatory_run_v01/seed_SEED.json, default 1)"""
import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
seed = int(sys.argv[1]) if len(sys.argv) > 1 else 1
with open(os.path.join(ROOT, "results", "confirmatory_run_v01", f"seed_{seed}.json")) as fh:
    rec = json.load(fh)

rows = [("rp_d0", "Radical-pair compass", "#534AB7"), ("polar", "Polar control", "#0F6E56")]
cols = [("Z2", "H1"), ("Z2", "H2"), ("Z3", "H1"), ("Z3", "H2")]
SLOTS = 8                                   # the seed files keep the 8 longest bars per dimension
fig, axes = plt.subplots(2, 4, figsize=(11, 5.0), sharex="col")
for i, (cond, label, colour) in enumerate(rows):
    r = rec["conditions"][cond]
    b = r["betti"]
    for j, (z, h) in enumerate(cols):
        ax = axes[i, j]
        bars = r["bars"][z][h]
        tau = r["thresholds"][z][h]
        ax.barh(list(range(len(bars))), bars, color=[colour if v > tau else "#B4B2A9" for v in bars], height=0.7)
        ax.axvline(tau, color="#5F5E5A", linestyle="--", linewidth=1)
        ax.set_ylim(SLOTS - 0.5, -0.5)
        ax.set_yticks([])
        if i == 0:
            ax.set_title(f"{h} with {z} coefficients", fontsize=11)
        if j == 0:
            ax.set_ylabel(f"{label}\nlabel: {r['shape']}", fontsize=10)
        if i == 1:
            ax.set_xlabel("bar lifetime", fontsize=10)
        ax.text(0.97, 0.06, f"Betti {b[z][h]}", transform=ax.transAxes, ha="right", va="bottom", fontsize=9,
                color="#444441")
        ax.spines[["top", "right"]].set_visible(False)
for j, (z, h) in enumerate(cols):                # room on the right of every threshold for the Betti note
    col_max = max(max(rec["conditions"][c]["bars"][z][h] + [rec["conditions"][c]["thresholds"][z][h]]) for c, _, _ in rows)
    axes[0, j].set_xlim(0, 1.3 * col_max)
fig.suptitle(f"Longest persistence bars, simulated bird {seed}; coloured bars pass the counting threshold (dashed)",
             fontsize=12)
fig.tight_layout()
out = os.path.join(ROOT, "figures", f"barcodes_seed{seed}.png")
fig.savefig(out, dpi=180, facecolor="white")
print(out)
