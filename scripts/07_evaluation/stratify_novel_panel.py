#!/usr/bin/env python3
"""
Stratify the ANI-novel 15-phage panel into "dark" vs "has_relatives" strata
and compute per-stratum recall for each tool at L1500 and L3000.

The 15-phage panel was selected using CheckV, which can only assign high quality
to phages with genes recognisable by its public-database-derived HMMs. This
creates a selection bias: genuinely dark phages are filtered out. Of the 15
selected phages, 13 retain protein-level relatives in MetaVR (84-95% pident),
while 2 (novel_Lcoc_phage_14, novel_Lcoc_phage_15) have no MetaVR hit at all.

This script stratifies the existing per-source-genome benchmark data to quantify
the impact of this selection bias on each tool's recall.

Input files (all pre-existing):
  - data/novel_spike_discovery/full_scale/novel_benchmark_L{1500,3000}_by_source_genome.tsv
  - data/novel_spike_discovery/full_scale/novel_genome_strata.tsv
    (short_id, metavr_best_pident; sanitised from the private selected/id_mapping.tsv)

Output:
  - data/novel_spike_discovery/full_scale/dark_vs_relatives_stratification.tsv

Usage:
    python stratify_novel_panel.py [--output PATH]
"""

import argparse
import csv
import os
import sys
from collections import defaultdict
from pathlib import Path


# ════════════════════════════════════════════════════════════════════
# Tool categories (methodological approaches; Table S21)
# ════════════════════════════════════════════════════════════════════

TOOL_CATEGORY = {
    "MetaPhinder": "A: Reference",
    "Sourmash": "A: Reference",
    "VirSorter": "B: ML+HMM",
    "VirFinder": "B: ML+HMM",
    "VirSorter2": "B: ML+HMM",
    "VIBRANT": "B: ML+HMM",
    "DeepVirFinder": "C: Deep Learning",
    "PPR-Meta": "C: Deep Learning",
    "Seeker": "C: Deep Learning",
    "Jaeger": "C: Deep Learning",
    "HVSeeker": "D: LM/Transformer",
    "TransGINmer": "D: LM/Transformer",
    "ViraLM": "D: LM/Transformer",
    "geNomad": "E: DNN+Markers",
}

CATEGORY_ORDER = [
    "A: Reference",
    "B: ML+HMM",
    "C: Deep Learning",
    "D: LM/Transformer",
    "E: DNN+Markers",
]


def load_stratum_map(id_mapping_path):
    """
    Load novel_genome_strata.tsv and classify each short_id as 'dark' or 'has_relatives'.

    A genome is 'dark' if metavr_best_pident == 'no_hit', meaning no nucleotide-
    level match to MetaVR v5. All other genomes are 'has_relatives'.
    """
    stratum = {}
    with open(id_mapping_path) as f:
        reader = csv.DictReader((l for l in f if not l.startswith("#")), delimiter="\t")
        for row in reader:
            sid = row["short_id"]
            pident = row["metavr_best_pident"]
            if pident == "no_hit":
                stratum[sid] = "dark"
            else:
                stratum[sid] = "has_relatives"
    return stratum


def load_by_source_genome(tsv_path, stratum_map):
    """
    Load per-source-genome benchmark TSV, filter to viral rows only,
    and annotate each row with its stratum.

    Returns list of dicts with added 'stratum' and 'length' keys.
    """
    rows = []
    with open(tsv_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            # Only viral rows (skip bacterial negative controls)
            if row["label"] != "viral":
                continue
            sg = row["source_genome"]
            if sg not in stratum_map:
                print(
                    f"WARNING: source_genome '{sg}' not found in novel_genome_strata.tsv, "
                    f"skipping",
                    file=sys.stderr,
                )
                continue
            row["stratum"] = stratum_map[sg]
            rows.append(row)
    return rows


def compute_stratified_metrics(rows, length_label):
    """
    Aggregate TP and FN per tool x stratum, compute recall.

    Returns list of result dicts.
    """
    # Accumulate: (tool, stratum) -> {TP, FN, n_fragments, n_genomes}
    agg = defaultdict(lambda: {"TP": 0, "FN": 0, "n_fragments": 0, "genomes": set()})

    for row in rows:
        tool = row["tool"]
        stratum = row["stratum"]
        key = (tool, stratum)
        tp = int(row["TP"])
        fn = int(row["FN"])
        agg[key]["TP"] += tp
        agg[key]["FN"] += fn
        agg[key]["n_fragments"] += int(row["n_fragments"])
        agg[key]["genomes"].add(row["source_genome"])

    results = []
    # Build per-tool recall by stratum
    tool_recalls = defaultdict(dict)  # tool -> {stratum: recall}

    for (tool, stratum), vals in sorted(agg.items()):
        tp = vals["TP"]
        fn = vals["FN"]
        total = tp + fn
        recall = tp / total if total > 0 else 0.0
        tool_recalls[tool][stratum] = recall
        results.append(
            {
                "category": TOOL_CATEGORY.get(tool, "Unknown"),
                "tool": tool,
                "stratum": stratum,
                "length": length_label,
                "n_genomes": len(vals["genomes"]),
                "n_fragments": vals["n_fragments"],
                "TP": tp,
                "FN": fn,
                "recall": round(recall, 4),
                "recall_delta_vs_relatives": None,  # filled below
            }
        )

    # Compute recall delta: dark - has_relatives
    for r in results:
        tool = r["tool"]
        if r["stratum"] == "dark" and "has_relatives" in tool_recalls.get(tool, {}):
            delta = tool_recalls[tool]["dark"] - tool_recalls[tool]["has_relatives"]
            r["recall_delta_vs_relatives"] = round(delta, 4)
        elif r["stratum"] == "has_relatives" and "dark" in tool_recalls.get(tool, {}):
            delta = tool_recalls[tool]["has_relatives"] - tool_recalls[tool]["dark"]
            r["recall_delta_vs_relatives"] = round(delta, 4)

    return results


def sort_results(results):
    """Sort by category order -> tool name -> stratum -> length."""
    cat_order = {c: i for i, c in enumerate(CATEGORY_ORDER)}
    return sorted(
        results,
        key=lambda r: (
            cat_order.get(r["category"], 99),
            r["tool"],
            r["stratum"],
            r["length"],
        ),
    )


def print_summary_table(results):
    """Print a formatted summary grouped by category."""
    print("\n" + "=" * 100)
    print("Dark vs Has-Relatives Stratification: Per-Tool Recall on ANI-Novel Panel")
    print("=" * 100)

    # Group by length
    for length in ["L1500", "L3000"]:
        print(f"\n{'─' * 100}")
        print(f"  {length}")
        print(f"{'─' * 100}")
        print(
            f"  {'Category':<20s} {'Tool':<16s} {'Stratum':<16s} "
            f"{'Genomes':>7s} {'Frags':>6s} {'TP':>5s} {'FN':>5s} "
            f"{'Recall':>7s} {'Delta':>8s}"
        )
        print(f"  {'─' * 95}")

        current_cat = None
        for r in results:
            if r["length"] != length:
                continue
            cat = r["category"]
            if cat != current_cat:
                if current_cat is not None:
                    print()
                current_cat = cat

            delta_str = (
                f"{r['recall_delta_vs_relatives']:+.4f}"
                if r["recall_delta_vs_relatives"] is not None
                else "  --"
            )
            print(
                f"  {cat:<20s} {r['tool']:<16s} {r['stratum']:<16s} "
                f"{r['n_genomes']:>7d} {r['n_fragments']:>6d} {r['TP']:>5d} "
                f"{r['FN']:>5d} {r['recall']:>7.4f} {delta_str:>8s}"
            )


def main():
    parser = argparse.ArgumentParser(
        description="Stratify ANI-novel panel into dark vs has_relatives strata"
    )
    parser.add_argument(
        "--base-dir",
        default="data/novel_spike_discovery/full_scale",
        help="Base directory containing novel benchmark TSVs and selected/ subdir",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Output TSV path (default: <base-dir>/dark_vs_relatives_stratification.tsv)",
    )
    args = parser.parse_args()

    base = Path(args.base_dir)
    id_mapping = base / "novel_genome_strata.tsv"
    l1500_tsv = base / "novel_benchmark_L1500_by_source_genome.tsv"
    l3000_tsv = base / "novel_benchmark_L3000_by_source_genome.tsv"

    # Validate inputs
    for p in [id_mapping, l1500_tsv, l3000_tsv]:
        if not p.exists():
            print(f"ERROR: required input file not found: {p}", file=sys.stderr)
            sys.exit(1)

    output_path = args.output or str(base / "dark_vs_relatives_stratification.tsv")

    # 1. Load stratum assignments
    stratum_map = load_stratum_map(id_mapping)
    n_dark = sum(1 for v in stratum_map.values() if v == "dark")
    n_rel = sum(1 for v in stratum_map.values() if v == "has_relatives")
    print(f"Loaded {len(stratum_map)} genomes: {n_dark} dark, {n_rel} has_relatives")

    # 2. Load per-source-genome TSVs (viral rows only)
    rows_l1500 = load_by_source_genome(l1500_tsv, stratum_map)
    rows_l3000 = load_by_source_genome(l3000_tsv, stratum_map)
    print(f"L1500: {len(rows_l1500)} viral rows loaded")
    print(f"L3000: {len(rows_l3000)} viral rows loaded")

    # 3. Compute stratified metrics
    results_l1500 = compute_stratified_metrics(rows_l1500, "L1500")
    results_l3000 = compute_stratified_metrics(rows_l3000, "L3000")

    all_results = sort_results(results_l1500 + results_l3000)

    # 4. Write output TSV
    fieldnames = [
        "category",
        "tool",
        "stratum",
        "length",
        "n_genomes",
        "n_fragments",
        "TP",
        "FN",
        "recall",
        "recall_delta_vs_relatives",
    ]
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for r in all_results:
            # Convert None to empty string for TSV
            row = dict(r)
            if row["recall_delta_vs_relatives"] is None:
                row["recall_delta_vs_relatives"] = ""
            writer.writerow(row)

    print(f"\nOutput written to: {output_path}")

    # 5. Print summary table
    print_summary_table(all_results)


if __name__ == "__main__":
    main()
