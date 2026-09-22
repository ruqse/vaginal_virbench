#!/usr/bin/env python3
"""
05_plot_vmgc_overlap.py — optional Part 2 figure.

Preferred visualization: a host-phylum x database ABSENCE view that keeps the
two databases' distinct semantics explicit (MetaVR = nucleotide species-level
presence; RefSeq = protein homology present). We deliberately AVOID a single
Venn that mixes unlike evidence.

Panel A: stacked bar of MetaVR status (species / relative-only / absent) for the
         top host phyla.
Panel B: % with NO RefSeq protein homolog, per host phylum (protein homology is a
         different axis -- labelled as such, not "catalogue overlap").

Input: results/tables/table_sXX_vmgc_benchmark_overlap.tsv (from 04).
"""

import argparse
import csv
import os
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# optional shared styling
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "07_evaluation"))
try:
    import pub_style
    pub_style.apply()
except Exception:  # noqa: BLE001
    pass

TOP_N_PHYLA = 8


def norm_phylum(p):
    """Collapse ambiguous multi-host predictions and unassigned hosts."""
    p = (p or "").strip()
    if p in ("", "-"):
        return "Unassigned"
    if "," in p:
        return "Mixed"
    return p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--overlap-table", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    by_phylum = defaultdict(lambda: {"present_metavr_species": 0,
                                     "present_metavr_relative_only": 0,
                                     "absent_metavr": 0,
                                     "refseq_homolog": 0, "n": 0})
    for_total = 0
    with open(a.overlap_table) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            for_total += 1
            p = norm_phylum(r["vmgc_host_phylum"])
            d = by_phylum[p]
            d["n"] += 1
            d[r["metavr_status"]] += 1
            if r["refseq_status"] == "refseq_protein_homolog":
                d["refseq_homolog"] += 1

    # Order: biological phyla by n (desc), then Mixed, then Unassigned last.
    special = {"Mixed", "Unassigned"}
    bio = sorted((p for p in by_phylum if p not in special),
                 key=lambda x: -by_phylum[x]["n"])[:TOP_N_PHYLA]
    phyla = bio + [p for p in ("Mixed", "Unassigned") if p in by_phylum]
    x = list(range(len(phyla)))

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(13, 5.2))

    species = [by_phylum[p]["present_metavr_species"] for p in phyla]
    relative = [by_phylum[p]["present_metavr_relative_only"] for p in phyla]
    absent = [by_phylum[p]["absent_metavr"] for p in phyla]
    ns = [by_phylum[p]["n"] for p in phyla]
    axA.bar(x, species, label="species match (ANI ≥ 95%)", color="#2c7fb8")
    axA.bar(x, relative, bottom=species, label="relative only (< species)", color="#7fcdbb")
    axA.bar(x, absent, bottom=[s + r for s, r in zip(species, relative)],
            label="no MetaVR v5 hit", color="#d95f0e")
    # annotate species-level ABSENCE (relative + absent) atop each bar
    for xi, p, n in zip(x, phyla, ns):
        no_species = by_phylum[p]["present_metavr_relative_only"] + by_phylum[p]["absent_metavr"]
        axA.text(xi, n + max(ns) * 0.01, f"{100.0*no_species/n:.0f}%",
                 ha="center", va="bottom", fontsize=8, color="#444444")
    axA.set_xticks(x)
    axA.set_xticklabels(phyla, rotation=35, ha="right")
    axA.set_ylabel("VMGC vOTUs")
    axA.set_ylim(0, max(ns) * 1.22)
    axA.set_title("A  MetaVR v5 (nucleotide, species-level)\n% = vOTUs lacking a species match")
    axA.legend(fontsize=8, frameon=False, loc="upper center", ncol=3,
               bbox_to_anchor=(0.5, 1.0), columnspacing=1.0, handletextpad=0.4)

    pct_refseq = [100.0 * by_phylum[p]["refseq_homolog"] / by_phylum[p]["n"] for p in phyla]
    axB.bar(x, pct_refseq, color="#54278f")
    for xi, v in zip(x, pct_refseq):
        axB.text(xi, v + 1, f"{v:.0f}", ha="center", va="bottom", fontsize=8, color="#444444")
    axB.set_xticks(x)
    axB.set_xticklabels(phyla, rotation=35, ha="right")
    axB.set_ylabel("% vOTUs with a RefSeq viral-protein homolog")
    axB.set_ylim(0, 105)
    axB.set_title("B  RefSeq viral protein (homology present)")

    fig.suptitle(
        f"VMGC vaginal vOTUs are nucleotide-novel to MetaVR v5 but retain protein homologs "
        f"(n = {for_total} vOTUs, by predicted host phylum)", fontsize=12)
    fig.tight_layout()
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    fig.savefig(a.out, dpi=300, bbox_inches="tight")
    print(f"  wrote {a.out}")


if __name__ == "__main__":
    main()
