#!/usr/bin/env python3
"""Break down which evidence lines contribute to each tier, and how many contigs.

For the >=1,500 bp benchmark window, for tiers 1 / 2 / 3:
  - per RAW evidence line (E1a,E1b,E2,E3,E4,E5): # contigs in the tier carrying it
  - GROUP-count distribution (the tier rule counts E1=E1a|E1b as one group): how
    many of the 5 groups {E1,E2,E3,E4,E5} fire
  - orthogonal (E4/E5) support
  - top evidence-line combinations
Outputs tables + a figure. Pure tabulation from the 13 per-sample GT files.

Run: module load SciPy-bundle/2024.05-gfbf-2024a matplotlib/3.9.2-gfbf-2024a
Out: results/.../analysis/gt_confidence/tier_evidence_*.tsv,
     results/figures/fig_s16.png, results/figures/fig_s16.pdf
"""
from __future__ import annotations
import csv
from pathlib import Path
from collections import Counter, defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
GTROOT = ROOT / "results/test_real/ground_truth"
OUT = ROOT / "results/test_real/viralm_investigation/analysis/gt_confidence"
FIG = ROOT / "results/figures"
OUT.mkdir(parents=True, exist_ok=True); FIG.mkdir(parents=True, exist_ok=True)
SAMPLES = ["UC028_V2", "UC055_V1", "UC055_V2", "UC062_V2", "UC065_V2", "UC074_V2",
           "UC084_V2", "UC093_V2", "UC093_V3", "UC096_V2", "UC115_V2", "UC139_V2", "UC164_V2"]
MIN_LEN = 1500
LINES = ["E1a", "E1b", "E2", "E3", "E4", "E5"]
ORTH = {"E4", "E5"}

# tier -> raw-line membership counts; tier -> group-count dist; tier -> combos
line_ct = defaultdict(Counter)
grp_dist = defaultdict(Counter)
combos = defaultdict(Counter)
n_tier = Counter()
orth_ct = Counter()
for s in SAMPLES:
    with open(GTROOT / s / "ground_truth_with_kraken2.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if int(r["length"]) < MIN_LEN:
                continue
            t = int(r["tier"])
            if t not in (1, 2, 3):
                continue
            ev = (r.get("evidence_lines", "-") or "-")
            L = set(ev.split(",")) if ev != "-" else set()
            n_tier[t] += 1
            grp_dist[t][int(r.get("n_evidence", 0))] += 1
            combos[t][ev] += 1
            if L & ORTH:
                orth_ct[t] += 1
            for x in LINES:
                if x in L:
                    line_ct[t][x] += 1

# ---- table 1: per-tier raw-line membership ----
rows = []
TIER_NAME = {1: "Tier 1 (>=3 groups, positive)", 2: "Tier 2 (=2 groups, positive)",
             3: "Tier 3 (=1 group, EXCLUDED)"}
for t in (1, 2, 3):
    n = n_tier[t]
    row = {"tier": TIER_NAME[t], "n_contigs": n}
    for x in LINES:
        c = line_ct[t][x]
        row[x] = f"{c} ({100*c/max(n,1):.0f}%)"
    row["orthogonal_E4|E5"] = f"{orth_ct[t]} ({100*orth_ct[t]/max(n,1):.0f}%)"
    rows.append(row)
with open(OUT / "tier_evidence_membership.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader(); w.writerows(rows)

# ---- table 2: group-count distribution per tier ----
with open(OUT / "tier_groupcount_distribution.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["tier", "n_groups (E1|E2|E3|E4|E5)", "n_contigs", "pct_of_tier"])
    for t in (1, 2, 3):
        for g in sorted(grp_dist[t]):
            w.writerow([TIER_NAME[t], g, grp_dist[t][g],
                        f"{100*grp_dist[t][g]/max(n_tier[t],1):.1f}"])

# ---- table 3: top combos per tier ----
with open(OUT / "tier_evidence_combos.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["tier", "evidence_lines", "n_contigs", "pct_of_tier"])
    for t in (1, 2, 3):
        for combo, c in combos[t].most_common():
            w.writerow([TIER_NAME[t], combo, c, f"{100*c/max(n_tier[t],1):.1f}"])

# ---- print summary ----
print("Per-tier evidence-line membership (>=1500 bp):")
for t in (1, 2, 3):
    print(f"  {TIER_NAME[t]}  n={n_tier[t]}")
    print("    " + "  ".join(f"{x}:{line_ct[t][x]}({100*line_ct[t][x]/max(n_tier[t],1):.0f}%)" for x in LINES)
          + f"   orth:{orth_ct[t]}({100*orth_ct[t]/max(n_tier[t],1):.0f}%)")
    print("    group-count: " + ", ".join(f"{g}grp:{grp_dist[t][g]}" for g in sorted(grp_dist[t])))

# ---- figure: heatmap (tier x line, % of tier) + group-count stacked bar ----
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"], "font.size": 9})
fig, (axH, axB) = plt.subplots(1, 2, figsize=(13, 4.2), gridspec_kw={"width_ratios": [1.55, 1.0], "wspace": 0.28})

tiers = [1, 2, 3]
M = np.array([[100*line_ct[t][x]/max(n_tier[t], 1) for x in LINES] for t in tiers])
Cnt = np.array([[line_ct[t][x] for x in LINES] for t in tiers])
im = axH.imshow(M, cmap="YlGnBu", vmin=0, vmax=100, aspect="auto")
axH.set_xticks(range(len(LINES))); axH.set_xticklabels(
    ["E1a\nprotein", "E1b\nMetaVR", "E2\nCheckV", "E3\nprophage", "E4\nCRISPR", "E5\nRCA"], fontsize=8.5)
axH.set_yticks(range(3))
axH.set_yticklabels([f"Tier 1\n(n={n_tier[1]})", f"Tier 2\n(n={n_tier[2]})", f"Tier 3\n(n={n_tier[3]})"], fontsize=9)
# group dividers: homology (E1a,E1b,E2,E3) vs orthogonal (E4,E5)
axH.axvline(3.5, color="#C44E52", lw=2)
axH.text(1.5, -0.75, "viral homology / genes", ha="center", color="#2c5a8c", fontsize=8.5, weight="bold")
axH.text(4.5, -0.75, "reference-free", ha="center", color="#C44E52", fontsize=8.5, weight="bold")
for i in range(3):
    for j in range(len(LINES)):
        axH.text(j, i, f"{Cnt[i,j]}\n{M[i,j]:.0f}%", ha="center", va="center",
                 fontsize=8, color=("white" if M[i, j] > 55 else "#222222"), weight="bold")
axH.set_title("Evidence-line membership by tier  (cell = # contigs, % of tier)", fontsize=10, weight="bold")
cb = fig.colorbar(im, ax=axH, fraction=0.046, pad=0.02); cb.set_label("% of tier's contigs", fontsize=8)

# group-count stacked bar: per tier, fraction by n_groups
maxg = 5
colors = {1: "#C44E52", 2: "#DD8452", 3: "#55A868", 4: "#4C72B0", 5: "#3a3a3a"}
bottoms = np.zeros(3)
for g in range(1, maxg + 1):
    vals = np.array([grp_dist[t].get(g, 0) for t in tiers], float)
    axB.bar(range(3), vals, bottom=bottoms, color=colors[g], edgecolor="white",
            label=f"{g} line" + ("s" if g > 1 else ""))
    for i, v in enumerate(vals):
        if v > 12:
            axB.text(i, bottoms[i] + v/2, f"{int(v)}", ha="center", va="center",
                     color="white", fontsize=8, weight="bold")
    bottoms += vals
axB.set_xticks(range(3)); axB.set_xticklabels([f"Tier 1\n(n={n_tier[1]})", f"Tier 2\n(n={n_tier[2]})", f"Tier 3\n(n={n_tier[3]})"])
axB.set_ylabel("contigs"); axB.spines[["top", "right"]].set_visible(False)
axB.set_title("Number of supporting evidence lines per tier\n(E1=E1a|E1b counted once)", fontsize=10, weight="bold")
axB.legend(fontsize=8, frameon=False, loc="upper left", title="lines: {E1,E2,E3,E4,E5}", title_fontsize=8)
fig.suptitle("Evidence-line contribution to each ground-truth tier (>= 1,500 bp)", fontsize=12.5, weight="bold", y=1.04)
# Canonical supplementary name only. The legacy "fig_tier_evidence_breakdown"
# PNG was byte-identical to fig_s16.png (retired 2026-08-11), and emitting the
# vector under the legacy name left Fig. S16 with no discoverable PDF.
for ext in ("png", "pdf"):
    fig.savefig(FIG / f"fig_s16.{ext}", dpi=300, bbox_inches="tight")
print("\nwrote", FIG / "fig_s16.png", "and", FIG / "fig_s16.pdf")
