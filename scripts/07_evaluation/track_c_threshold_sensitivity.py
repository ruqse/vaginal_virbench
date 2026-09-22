#!/usr/bin/env python3
"""
track_c_threshold_sensitivity.py — Enrichment-ratio (R) sensitivity for Track C

Complements the RCA mapping-support grid in rca_sensitivity_grid.py (which
sweeps breadth × read_pairs) by sweeping the *R* enrichment-ratio thresholds
that define Track C TP / TN / dark-matter labels.

Current Track C defaults (from compute_enrichment.sh and build_rca_ground_truth.py):
    TP: R > 10 AND CheckV viral signal
    TN: 0.5 < R < 2 AND Kraken2 bacterial AND CheckV viral_genes = 0
    dark matter: R > 10 AND no CheckV signal AND no host genes AND unclassified
    excluded: everything else

This script re-labels contigs at a grid of alternative R cutoffs without
rerunning CheckV, Kraken2, or the 14 tools — it reuses:
  - *_enrichment.tsv (raw R values per contig, 13 patients)
  - checkv/quality_summary.tsv (pooled, single file)
  - kraken2_output.tsv + kraken2_report.txt
  - per-tool predictions already parsed for the default Track C benchmark

At each R combination it computes per-tool MCC on the resulting TP / TN set
and writes:
  - track_c_sensitivity_mcc.tsv   (one row per threshold combo, tool columns)
  - track_c_sensitivity_counts.tsv (TP/TN counts per combo)
  - track_c_sensitivity_ranking.tsv (top-5 tools per combo)

Usage:
    python track_c_threshold_sensitivity.py \\
        --coassembly-dir results/test_real/coassembly \\
        --checkv-dir     results/test_real/coassembly/checkv \\
        --kraken2-output results/test_real/coassembly/kraken2/kraken2_output.tsv \\
        --results-dir    results/test_real/full_run/RCA_master \\
        --output-dir     results/test_real/coassembly/sensitivity_R
"""

import argparse
import csv
import importlib.util
import sys
from collections import Counter
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Import parsers + helpers from existing modules (don't reimplement)
# ---------------------------------------------------------------------------

_SCRIPTS = Path(__file__).resolve().parent.parent
_BUILD_GT = _SCRIPTS / "03_assembly" / "build_rca_ground_truth.py"
_EVAL_RCA = _SCRIPTS / "03_assembly" / "evaluate_rca_benchmark.py"


def _load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_bgt = _load_module("build_rca_ground_truth", _BUILD_GT)
_eval = _load_module("evaluate_rca_benchmark", _EVAL_RCA)

PATIENTS = _bgt.PATIENTS
BACTERIAL_DOMAINS = _bgt.BACTERIAL_DOMAINS
has_viral_signal = _bgt.has_viral_signal
load_enrichment = _bgt.load_enrichment
load_checkv = _bgt.load_checkv
load_kraken2 = _bgt.load_kraken2
classify_domain = _bgt.classify_domain

TOOL_THRESHOLDS = _eval.TOOL_THRESHOLDS
TOOL_PVALUE_THRESHOLDS = _eval.TOOL_PVALUE_THRESHOLDS
compute_metrics = _eval.compute_metrics
evaluate_tool = _eval.evaluate_tool

TOOLS = list(TOOL_THRESHOLDS.keys())


# ---------------------------------------------------------------------------
# Threshold grid
# ---------------------------------------------------------------------------

# TP threshold: R > r_enriched_min
TP_THRESHOLDS = [5.0, 10.0, 20.0]

# TN band: bg_low < R < bg_high
#   (0.33, 3.00) wide symmetric band around R = 1
#   (0.50, 2.00) current default
#   (0.80, 1.25) tight symmetric band
TN_BANDS = [(0.33, 3.00), (0.50, 2.00), (0.80, 1.25)]

DEFAULT_TP = 10.0
DEFAULT_TN = (0.50, 2.00)


# ---------------------------------------------------------------------------
# Per-combination ground truth builder (in-memory, no disk I/O per combo)
# ---------------------------------------------------------------------------

def build_gt_at_threshold(enrichment, checkv, kraken, domains,
                          r_tp_min, r_bg_low, r_bg_high):
    """Re-classify contigs at the given R thresholds.

    Returns {contig_id: {'label': 'viral'|'negative', ...}}, TP/TN/dark counts.
    Only TP and TN are returned; dark_matter and excluded are dropped.
    """
    gt = {}
    counts = Counter()

    for contig_id in sorted(enrichment.keys()):
        e = enrichment[contig_id]
        r = e["enrichment_R"]
        cv = checkv.get(contig_id, {
            "checkv_quality": "Not-determined",
            "viral_genes": 0,
            "host_genes": 0,
            "provirus": "No",
        })
        kr = kraken.get(contig_id, {"classified": False, "taxon_name": "unclassified"})
        domain = domains.get(contig_id, "unknown")

        is_viral, _ = has_viral_signal(cv)

        # TP: R > r_tp_min AND CheckV viral signal
        if r > r_tp_min:
            if is_viral:
                gt[contig_id] = {
                    "label": "viral",
                    "contig_length": e["contig_length"],
                    "checkv_quality": cv.get("checkv_quality", "NA"),
                    "viral_genes": cv.get("viral_genes", 0),
                    "host_genes": cv.get("host_genes", 0),
                    "provirus": cv.get("provirus", "No"),
                }
                counts["TP"] += 1
            elif (cv.get("viral_genes", 0) == 0
                  and cv.get("host_genes", 0) == 0
                  and not kr["classified"]):
                counts["dark_matter"] += 1
            else:
                counts["excluded"] += 1
        # TN: bg_low < R < bg_high AND bacterial AND no viral signal
        elif r_bg_low < r < r_bg_high:
            if not is_viral and domain in BACTERIAL_DOMAINS:
                gt[contig_id] = {
                    "label": "negative",
                    "contig_length": e["contig_length"],
                    "checkv_quality": cv.get("checkv_quality", "NA"),
                    "viral_genes": cv.get("viral_genes", 0),
                    "host_genes": cv.get("host_genes", 0),
                    "provirus": cv.get("provirus", "No"),
                }
                counts["TN"] += 1
            else:
                counts["excluded"] += 1
        else:
            counts["excluded"] += 1

    return gt, counts


# ---------------------------------------------------------------------------
# Tool prediction loader
# ---------------------------------------------------------------------------

def load_all_tool_predictions(results_dir: Path, sample_id: str = "RCA_master",
                              gt_ids=None) -> dict:
    """Load raw scores for all 14 tools once. Returns {tool: {contig: score_or_(score,p)}}."""
    scores = _eval.discover_tool_outputs(results_dir, sample_id, gt_ids=gt_ids)
    return {tool: scores.get(tool, {}) for tool in TOOLS}


# ---------------------------------------------------------------------------
# Main sweep
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--coassembly-dir", required=True, type=Path)
    parser.add_argument("--checkv-dir", required=True, type=Path)
    parser.add_argument("--kraken2-output", required=True, type=Path)
    parser.add_argument("--kraken2-report", type=Path, default=None)
    parser.add_argument("--results-dir", required=True, type=Path,
                        help="Root dir of per-tool predictions for RCA_master sample")
    parser.add_argument("--sample-id", default="RCA_master")
    parser.add_argument("--output-dir", "-o", required=True, type=Path)
    parser.add_argument("--tp-thresholds", nargs="+", type=float, default=TP_THRESHOLDS)
    parser.add_argument("--tn-bands", nargs="+", default=None,
                        help='Optional "low:high" pairs, e.g. "0.33:3.0 0.5:2.0"')
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    tp_vals = list(args.tp_thresholds)
    tn_bands = TN_BANDS
    if args.tn_bands:
        tn_bands = [tuple(float(x) for x in b.split(":")) for b in args.tn_bands]

    print("=" * 70)
    print("Track C enrichment-ratio threshold sensitivity")
    print("=" * 70)
    print(f"TP thresholds: {tp_vals}")
    print(f"TN bands:      {tn_bands}")
    print(f"Grid size:     {len(tp_vals) * len(tn_bands)} combinations")

    # --- Load data once -----------------------------------------------------
    print("\nLoading enrichment / CheckV / Kraken2...")
    enrichment = load_enrichment(args.coassembly_dir)
    checkv = load_checkv(args.checkv_dir)
    kraken = load_kraken2(args.kraken2_output)
    domains = classify_domain(args.kraken2_output, args.kraken2_report)
    print(f"  enrichment: {len(enrichment)}, checkv: {len(checkv)}, "
          f"kraken2: {len(kraken)}, domains: {len(domains)}")

    # --- Load all tool predictions once -------------------------------------
    print(f"\nLoading tool predictions from {args.results_dir}...")
    tool_preds = load_all_tool_predictions(args.results_dir, args.sample_id,
                                           gt_ids=set(enrichment))

    # --- Sweep grid ---------------------------------------------------------
    print("\nSweeping threshold grid...")
    mcc_rows = []
    count_rows = []
    ranking_rows = []

    for r_tp in tp_vals:
        for r_lo, r_hi in tn_bands:
            gt, counts = build_gt_at_threshold(
                enrichment, checkv, kraken, domains,
                r_tp_min=r_tp, r_bg_low=r_lo, r_bg_high=r_hi,
            )
            n_tp = counts["TP"]
            n_tn = counts["TN"]
            print(f"  R>{r_tp}, TN=({r_lo},{r_hi}):  "
                  f"TP={n_tp}, TN={n_tn}, dark={counts['dark_matter']}, "
                  f"excl={counts['excluded']}")

            is_default = (r_tp == DEFAULT_TP
                          and (r_lo, r_hi) == DEFAULT_TN)

            if n_tp == 0 or n_tn == 0:
                # Cannot compute MCC with no positives or no negatives
                mcc_rows.append({
                    "r_tp_min": r_tp, "r_bg_low": r_lo, "r_bg_high": r_hi,
                    "n_TP": n_tp, "n_TN": n_tn, "is_default": int(is_default),
                    **{t: "NA" for t in TOOLS},
                })
                count_rows.append({
                    "r_tp_min": r_tp, "r_bg_low": r_lo, "r_bg_high": r_hi,
                    "n_TP": n_tp, "n_TN": n_tn,
                    "n_dark_matter": counts["dark_matter"],
                    "n_excluded": counts["excluded"],
                    "is_default": int(is_default),
                })
                continue

            # Score each tool
            per_tool = {}
            for tool in TOOLS:
                threshold = TOOL_THRESHOLDS[tool]
                metrics = evaluate_tool(tool, tool_preds[tool], gt, threshold)
                per_tool[tool] = metrics

            mcc_row = {
                "r_tp_min": r_tp, "r_bg_low": r_lo, "r_bg_high": r_hi,
                "n_TP": n_tp, "n_TN": n_tn, "is_default": int(is_default),
            }
            for tool in TOOLS:
                mcc_row[tool] = per_tool[tool]["MCC"]
            mcc_rows.append(mcc_row)

            count_rows.append({
                "r_tp_min": r_tp, "r_bg_low": r_lo, "r_bg_high": r_hi,
                "n_TP": n_tp, "n_TN": n_tn,
                "n_dark_matter": counts["dark_matter"],
                "n_excluded": counts["excluded"],
                "is_default": int(is_default),
            })

            ranked = sorted(TOOLS, key=lambda t: -per_tool[t]["MCC"])
            ranking = {
                "r_tp_min": r_tp, "r_bg_low": r_lo, "r_bg_high": r_hi,
                "n_TP": n_tp, "n_TN": n_tn, "is_default": int(is_default),
            }
            for i, t in enumerate(ranked[:5]):
                ranking[f"rank_{i+1}"] = f"{t} ({per_tool[t]['MCC']:.4f})"
            ranking_rows.append(ranking)

    # --- Write outputs ------------------------------------------------------
    def write(rows, name, fields):
        p = args.output_dir / name
        with open(p, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields, delimiter="\t",
                               extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)
        print(f"  wrote {p}")

    mcc_fields = ["r_tp_min", "r_bg_low", "r_bg_high",
                  "n_TP", "n_TN", "is_default"] + TOOLS
    write(mcc_rows, "track_c_sensitivity_mcc.tsv", mcc_fields)

    cnt_fields = ["r_tp_min", "r_bg_low", "r_bg_high",
                  "n_TP", "n_TN", "n_dark_matter", "n_excluded", "is_default"]
    write(count_rows, "track_c_sensitivity_counts.tsv", cnt_fields)

    rank_fields = ["r_tp_min", "r_bg_low", "r_bg_high", "n_TP", "n_TN",
                   "is_default"] + [f"rank_{i+1}" for i in range(5)]
    write(ranking_rows, "track_c_sensitivity_ranking.tsv", rank_fields)

    # --- Ranking stability summary -----------------------------------------
    print("\n" + "=" * 70)
    print("Ranking stability")
    print("=" * 70)
    top3_orderings = [tuple(r[f"rank_{i+1}"].split(" (")[0] for i in range(3))
                      for r in ranking_rows]
    top3_sets = [frozenset(t) for t in top3_orderings]
    n = len(ranking_rows)
    c_ord = Counter(top3_orderings).most_common(1)[0]
    c_set = Counter(top3_sets).most_common(1)[0]
    print(f"  combos scored: {n}")
    print(f"  most common top-3 ordering: {c_ord[0]}  ({c_ord[1]}/{n})")
    print(f"  most common top-3 set (any order): "
          f"{sorted(c_set[0])}  ({c_set[1]}/{n})")


if __name__ == "__main__":
    main()
