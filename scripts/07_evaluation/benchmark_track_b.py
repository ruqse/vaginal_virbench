#!/usr/bin/env python3
"""
Benchmark Track B (assembly-based) virus identification tool predictions.

Evaluates 14 tools on Track B co-assemblies using per-assembly ground truth
(viral/bacterial/background labels from label_spike_in_contigs.py).

Ground truth interpretation (per discussion on Track B design):
  - "viral" contigs (mapped to spike-in genomes) = POSITIVE
  - "bacterial" + "background" contigs = NEGATIVE
  - "excluded_chimeric_risk" contigs = EXCLUDED (vB_Gva_AB1 on UC093_V3)

This differs from Track A where bacterial negatives are from fragmented reference
genomes. In Track B, the real metagenome itself serves as the negative set.

IMPORTANT — Non-monotonic evaluation set size:
    The number of viral contigs (n_viral) is non-monotonic with coverage depth
    due to metaSPAdes assembly dynamics:
      - 0.1x-1x: 0 viral contigs (below assembly threshold)
      - 5x: marginal assembly (UC115_V2: 44, UC093_V3: 35), mostly <2 kb
      - 10x: peak diversity (UC115_V2: 117, UC093_V3: 103), all categories
      - 50x: assembly collapse (both: 31), few long herpesvirus-dominated contigs
    Metrics at different coverages are therefore computed on evaluation sets that
    differ in size AND composition.  Cross-coverage comparisons are confounded
    by this changing composition.  The 10x condition is the most informative
    benchmark point.
    See: scripts/07_evaluation/plot_track_b_assembly_dynamics.py

Outputs:
  - Per-assembly metrics TSV
  - Aggregated summary across coverages
  - Recall-by-coverage curves

Usage:
    python benchmark_track_b.py \
        --assemblies-dir data/spike_in/assemblies/ \
        --results-dir results/spike_in_benchmark/ \
        --output-dir results/spike_in_benchmark/trackB_metrics/
"""

import argparse
import csv
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

# ════════════════════════════════════════════════════════════════════
# Import tool parsers from benchmark_fragments.py (same tool output formats)
# ════════════════════════════════════════════════════════════════════

# Add parent dir to path to import from benchmark_fragments
sys.path.insert(0, str(Path(__file__).parent))
from benchmark_fragments import (
    TOOL_THRESHOLDS as FRAGMENT_TOOL_THRESHOLDS,
    TOOL_PVALUE_THRESHOLDS,
    parse_deepvirfinder, parse_virsorter, parse_virsorter2, parse_virfinder,
    parse_pprmeta, parse_vibrant, parse_seeker, parse_metaphinder,
    parse_sourmash, parse_genomad, parse_hvseeker, parse_jaeger,
    parse_transginmer, parse_viralm,
    compute_metrics,
)

# Retain the original operating cutoffs for comparability. Use an independent
# copy so Track B configuration cannot mutate another benchmark track.
TOOL_THRESHOLDS = {
    **FRAGMENT_TOOL_THRESHOLDS,
    "Sourmash": 0.5,
    # The shared parser maps categories 1, 2, 3, 4 to 1.0, .85, .70, .55.
    "VirSorter": 0.5,
}

TOOL_PARSERS = {
    "DeepVirFinder": parse_deepvirfinder,
    "VirSorter": parse_virsorter,
    "VirSorter2": parse_virsorter2,
    "VirFinder": parse_virfinder,
    "PPR-Meta": parse_pprmeta,
    "VIBRANT": parse_vibrant,
    "Seeker": parse_seeker,
    "MetaPhinder": parse_metaphinder,
    "Sourmash": parse_sourmash,
    "geNomad": parse_genomad,
    "HVSeeker": parse_hvseeker,
    "Jaeger": parse_jaeger,
    "TransGINmer": parse_transginmer,
    "ViraLM": parse_viralm,
}


def load_ground_truth(gt_path):
    """
    Load Track B ground truth.
    Returns dict: {contig_id: {"label": "viral"/"negative"/"excluded", "category": str}}
    """
    gt = {}
    with open(gt_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            cid = row["contig_id"]
            label = row["label"]
            category = row.get("source_category", "-")

            if label == "viral":
                gt[cid] = {"label": "viral", "category": category}
            elif label == "excluded_chimeric_risk":
                gt[cid] = {"label": "excluded", "category": category}
            else:
                # bacterial + background → negative
                gt[cid] = {"label": "negative", "category": category}
    return gt


def resolve_virsorter_assembly_ids(predictions, ground_truth_ids):
    """Recover assembly IDs after VirSorter replaces dots with underscores.

    Match complete IDs, preserving already-correct IDs and leaving unrelated
    predictions unmatched. Ambiguous normalization must not silently assign a
    call to the wrong contig.
    """
    reverse_ids = {}
    for cid in ground_truth_ids:
        mangled = cid.replace(".", "_")
        if mangled in reverse_ids and reverse_ids[mangled] != cid:
            raise ValueError(
                f"VirSorter ID normalization collision: {reverse_ids[mangled]!r} "
                f"and {cid!r} both map to {mangled!r}"
            )
        reverse_ids[mangled] = cid

    def exact_id(candidate):
        return candidate if candidate in ground_truth_ids else reverse_ids.get(candidate)

    resolved = {}
    for cid, score in predictions.items():
        original = exact_id(cid)
        # Strip only recognized VirSorter annotations, and only when the
        # complete identifier was not already recognized as a true contig ID.
        if original is None:
            body = cid.removesuffix("-circular")
            original = exact_id(body)
            if original is None:
                region = re.fullmatch(r"(.+)_gene_\d+_gene_\d+-\d+-\d+", body)
                if region:
                    parent = region.group(1)
                    original = exact_id(parent) or exact_id(parent.removesuffix("-circular"))
        if original is None:
            original = cid
        # Multiple regions may derive from one parent contig. Its strongest
        # category determines the call and its ranking score.
        resolved[original] = max(score, resolved.get(original, score))
    return resolved


def evaluate_assembly(gt, tool_results_dir, tool_parsers, tool_thresholds):
    """
    Evaluate all tools on one assembly.
    Returns list of per-tool result dicts.
    """
    results = []
    viral_ids = {cid for cid, v in gt.items() if v["label"] == "viral"}
    negative_ids = {cid for cid, v in gt.items() if v["label"] == "negative"}
    excluded_ids = {cid for cid, v in gt.items() if v["label"] == "excluded"}

    # Benchmarkable set = viral + negative (exclude "excluded")
    bench_ids = sorted(viral_ids | negative_ids)
    y_true = np.array([1 if cid in viral_ids else 0 for cid in bench_ids])

    for tool_name, parser in tool_parsers.items():
        threshold = tool_thresholds.get(tool_name, 0.5)

        try:
            if tool_name in {"VirSorter", "VirSorter2", "VIBRANT", "geNomad"}:
                predictions = parser(Path(tool_results_dir), manifest_ids=set(gt))
            else:
                predictions = parser(Path(tool_results_dir))
            if tool_name == "VirSorter":
                predictions = resolve_virsorter_assembly_ids(predictions, set(gt))
        except Exception as exc:
            raise RuntimeError(
                f"{tool_name}: failed to parse predictions in {tool_results_dir}"
            ) from exc

        # AUPRC uses the viral score alone. A p-value is a separate acceptance
        # criterion, never a viral score; larger p-values cannot create calls.
        present = np.array([cid in predictions for cid in bench_ids], dtype=bool)
        pvalue_threshold = TOOL_PVALUE_THRESHOLDS.get(tool_name)
        if pvalue_threshold is not None:
            raw_scores = [predictions.get(cid, (0.0, 1.0)) for cid in bench_ids]
            y_scores = np.array([float(value[0]) for value in raw_scores])
            y_pvalues = np.array([float(value[1]) for value in raw_scores])
            y_pred = (present & (y_scores >= threshold)
                      & (y_pvalues < pvalue_threshold)).astype(int)
        else:
            raw_scores = [predictions.get(cid, 0.0) for cid in bench_ids]
            y_scores = np.array([
                float(value) if np.isscalar(value) else float(np.max(value))
                for value in raw_scores
            ])
            y_pred = (present & (y_scores >= threshold)).astype(int)

        metrics = compute_metrics(y_true, y_pred, y_scores)
        # A nonfinite reported score does not define a valid ranking. Preserve
        # the existing binary calls, but do not infer a replacement AP score.
        if not np.all(np.isfinite(y_scores)):
            metrics["AUPRC"] = "NA"

        # Per-category recall (viral categories only)
        cat_counts = defaultdict(lambda: {"tp": 0, "fn": 0})
        for cid, positive in zip(bench_ids, y_pred):
            if cid not in viral_ids:
                continue
            cat = gt[cid]["category"]
            if positive:
                cat_counts[cat]["tp"] += 1
            else:
                cat_counts[cat]["fn"] += 1

        cat_recall = {}
        for cat, counts in cat_counts.items():
            total = counts["tp"] + counts["fn"]
            cat_recall[cat] = counts["tp"] / total if total > 0 else 0.0

        results.append({
            "tool": tool_name,
            "n_total": len(bench_ids),
            "n_viral": len(viral_ids),
            "n_negative": len(negative_ids),
            "n_excluded": len(excluded_ids),
            "n_predictions": len(predictions),
            **metrics,
            "category_recall": cat_recall,
        })

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Benchmark Track B assembly-based tool predictions.")
    parser.add_argument("--assemblies-dir", required=True, type=Path,
                        help="Path to spike_in/assemblies/ (contains UC115_V2/, UC093_V3/, control/)")
    parser.add_argument("--results-dir", required=True, type=Path,
                        help="Path to results/spike_in_benchmark/ (contains trackB_* dirs)")
    parser.add_argument("--output-dir", required=True, type=Path,
                        help="Output directory for Track B metrics")
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    backgrounds = ["UC115_V2", "UC093_V3"]
    coverages = ["0.1", "0.5", "1", "5", "10", "50"]

    all_results = []

    for bg in backgrounds:
        print(f"\n{'='*60}")
        print(f"Background: {bg}")
        print(f"{'='*60}")

        for cov in coverages:
            # Ground truth
            gt_path = args.assemblies_dir / bg / f"assembly_cov{cov}" / "ground_truth.tsv"
            if not gt_path.exists():
                print(f"  cov{cov}: SKIP (no ground truth)")
                continue

            gt = load_ground_truth(gt_path)
            n_viral = sum(1 for v in gt.values() if v["label"] == "viral")
            n_neg = sum(1 for v in gt.values() if v["label"] == "negative")
            n_excl = sum(1 for v in gt.values() if v["label"] == "excluded")

            if n_viral == 0:
                print(f"  cov{cov}: SKIP (0 viral contigs)")
                continue

            # Tool results
            tool_dir = args.results_dir / f"trackB_{bg}_cov{cov}"
            if not tool_dir.exists():
                print(f"  cov{cov}: SKIP (no tool results)")
                continue

            print(f"  cov{cov}: {n_viral} viral, {n_neg} negative, {n_excl} excluded")

            results = evaluate_assembly(gt, tool_dir, TOOL_PARSERS, TOOL_THRESHOLDS)

            for r in results:
                r["background"] = bg
                r["coverage"] = cov
                all_results.append(r)

    if not all_results:
        print("\nERROR: No results computed. Check paths.")
        sys.exit(1)

    # ── Write per-assembly metrics ──────────────────────────────────
    out_path = args.output_dir / "trackB_metrics.tsv"
    fieldnames = ["background", "coverage", "tool", "n_total", "n_viral", "n_negative",
                  "n_excluded", "n_predictions", "TP", "FP", "TN", "FN",
                  "precision", "recall", "F1", "MCC", "AUPRC"]

    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        for r in all_results:
            writer.writerow(r)

    print(f"\nMetrics written: {out_path}")

    # ── Summary table: recall by tool × coverage ───────────────────
    print(f"\n{'='*80}")
    print(f"Recall by Tool × Coverage")
    print(f"{'='*80}")

    for bg in backgrounds:
        print(f"\n--- {bg} ---")
        bg_results = [r for r in all_results if r["background"] == bg]
        tools = sorted(set(r["tool"] for r in bg_results))
        covs = sorted(set(r["coverage"] for r in bg_results), key=lambda x: float(x))

        header = f"{'Tool':<18}" + "".join(f"{'cov'+c:>10}" for c in covs)
        print(header)
        print("-" * len(header))

        for tool in tools:
            vals = []
            for cov in covs:
                match = [r for r in bg_results if r["tool"] == tool and r["coverage"] == cov]
                if match:
                    vals.append(f"{match[0]['recall']:>10.3f}")
                else:
                    vals.append(f"{'N/A':>10}")
            print(f"{tool:<18}" + "".join(vals))

    # ── Write recall-by-coverage for plotting ──────────────────────
    recall_path = args.output_dir / "trackB_recall_by_coverage.tsv"
    with open(recall_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["background", "coverage", "tool", "recall",
                                                "precision", "MCC", "F1", "n_viral"],
                                delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for r in all_results:
            writer.writerow(r)

    print(f"\nRecall-by-coverage: {recall_path}")
    print(f"\nDone. Generate figures with plot_track_b.py")


if __name__ == "__main__":
    main()
