#!/usr/bin/env python3
"""
11_tool_calls_on_plasmid_like.py — do the 14 tools call plasmid-like sequence viral?

Compares per-tool positive-call rates (default thresholds, Track C rules) on:
  DM        : the 48 dark-matter contigs
  TP        : Track C true positives (n = 29)
  TN_plas   : Track C true negatives that geNomad reports as plasmids
              (Kraken2 Bacteria/Archaea, no CheckV viral genes, 0.5 < R < 2)
  TN_chrom  : remaining Track C true negatives
  BACT_plas : any contig (not DM/TP) with Kraken2 Bacteria, 0 CheckV viral
              genes, and a geNomad plasmid call, regardless of R

Caveat: "plasmid" here is geNomad's classifier call, not an independent
plasmid-database match.
"""
import csv
import os
import sys
from collections import defaultdict
from pathlib import Path

PROJ = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJ / "scripts/03_assembly"))
import evaluate_rca_benchmark as ev  # noqa: E402

GT = PROJ / "results/test_real/coassembly/rca_ground_truth.tsv"
RES = PROJ / "results/test_real/full_run/RCA_master"
GEN = RES / "genomad/RCA_master_contigs_summary/RCA_master_contigs_plasmid_summary.tsv"
OUT = PROJ / "results/vmgc_crosscheck/dark_matter_metavr_refseq/tool_calls_by_plasmid_group.tsv"

rows = {r["contig_id"]: r for r in csv.DictReader(open(GT), delimiter="\t")}
plas = {r["seq_name"] for r in csv.DictReader(open(GEN), delimiter="\t")}

groups = defaultdict(set)
for cid, r in rows.items():
    lab = r["label"]
    if lab == "dark_matter":
        groups["DM"].add(cid)
    elif lab == "TP":
        groups["TP"].add(cid)
    elif lab == "TN":
        groups["TN_plas" if cid in plas else "TN_chrom"].add(cid)
    if (lab not in ("dark_matter", "TP") and cid in plas
            and r["kraken2_domain"] == "Bacteria" and r["viral_genes"] == "0"):
        groups["BACT_plas"].add(cid)

all_ids = set(rows)
scores = ev.discover_tool_outputs(RES, "RCA_master", gt_ids=all_ids)

order = ["DM", "TP", "TN_plas", "TN_chrom", "BACT_plas"]
out = []
for tool in sorted(scores):
    s = ev.fold_scores(scores[tool], tool, all_ids)
    thr = ev.TOOL_THRESHOLDS.get(tool, 0.5)
    pthr = ev.TOOL_PVALUE_THRESHOLDS.get(tool)

    def pos(c):
        v = s.get(c)
        if v is None:
            return 0
        if isinstance(v, tuple):
            return int(v[0] >= thr and (pthr is None or v[1] < pthr))
        return int(float(v) >= thr)

    rec = {"tool": tool}
    for g in order:
        n = len(groups[g])
        k = sum(pos(c) for c in groups[g])
        rec[g] = f"{k}/{n} ({100 * k / n:.0f}%)" if n else "NA"
    out.append(rec)

with open(OUT, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["tool"] + order, delimiter="\t")
    w.writeheader()
    w.writerows(out)
print({g: len(groups[g]) for g in order})
print(f"wrote {OUT}")
