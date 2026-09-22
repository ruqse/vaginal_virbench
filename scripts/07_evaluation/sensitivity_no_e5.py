#!/usr/bin/env python3
"""Per-sample evaluation of all 14 tools against the primary GT and the
no_E5 GT variant, pooled across the 13 matched-virome samples.

Emits a tool-level delta table (MCC and rank under each GT) suitable for
Supplementary Table S2 as the third sensitivity row.

Usage:
    python scripts/07_evaluation/sensitivity_no_e5.py \\
        --samples UC084_V2 UC115_V2 UC164_V2 UC093_V2 UC139_V2 \\
                  UC055_V1 UC065_V2 UC096_V2 UC028_V2 \\
                  UC074_V2 UC093_V3 UC055_V2 UC062_V2 \\
        --results-dir results/test_real/full_run \\
        --gt-root     results/test_real/ground_truth \\
        --output-dir  results/test_real/sensitivity_no_e5
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import pandas as pd

# Reuse parsers and thresholds from evaluate_metagenome.py
sys.path.insert(0, str(Path(__file__).resolve().parent))
from evaluate_metagenome import (  # noqa: E402
    TOOL_PVALUE_THRESHOLDS,
    TOOL_THRESHOLDS,
    discover_tool_outputs,
    load_ground_truth,
)

import csv
import numpy as np


def _parse_jaeger_recursive(jaeger_dir: Path) -> dict:
    """Recursive Jaeger parser — outputs live at jaeger/<sample>_contigs/*.tsv
    which the upstream non-recursive glob in evaluate_metagenome.py misses."""
    scores: dict = {}
    if not jaeger_dir.exists():
        return scores
    for tsv in jaeger_dir.rglob("*_default_jaeger.tsv"):
        with open(tsv) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                contig = row.get("contig_id", "").strip()
                # Jaeger rewrites the merged sample-prefixed id separator "__"
                # into a comma (e.g. MITCH01__NODE_7 -> MITCH01,NODE_7). Restore
                # it so keys match the GT. No-op for un-prefixed (comma-free) ids.
                if "," in contig:
                    contig = contig.replace(",", "__")
                try:
                    logits = np.array([
                        float(row["bacteria_score"]),
                        float(row["phage_score"]),
                        float(row["eukarya_score"]),
                        float(row["archaea_score"]),
                    ])
                except (KeyError, ValueError):
                    continue
                exp_l = np.exp(logits - logits.max())
                p_phage = float(exp_l[1] / exp_l.sum())
                if contig:
                    scores[contig] = p_phage
    return scores


def _predict(score, threshold, pvalue_thr):
    """Return True iff the tool predicts 'viral' under (threshold, pvalue_thr).

    Mirrors evaluate_tool() in evaluate_metagenome.py: score values may be bare
    floats, (score, pvalue) tuples, or missing (→ 0.0).
    """
    if isinstance(score, tuple):
        s, pv = score[0], score[1]
        if pvalue_thr is not None:
            return (s >= threshold) and (pv < pvalue_thr)
        return s >= threshold
    return float(score) >= threshold


def confusion(tool_scores: dict, gt: dict, threshold: float,
              pvalue_thr) -> tuple[int, int, int, int]:
    """TP/FP/TN/FN over gt contigs. Contigs absent from tool_scores → negative."""
    tp = fp = tn = fn = 0
    for cid, info in gt.items():
        label = info["label"]
        score = tool_scores.get(cid, 0.0)
        pred_pos = _predict(score, threshold, pvalue_thr)
        if label == "viral":
            if pred_pos:
                tp += 1
            else:
                fn += 1
        elif label == "negative":
            if pred_pos:
                fp += 1
            else:
                tn += 1
    return tp, fp, tn, fn


def mcc(tp: int, fp: int, tn: int, fn: int) -> float:
    denom = math.sqrt(
        max((tp + fp), 1)
        * max((tp + fn), 1)
        * max((tn + fp), 1)
        * max((tn + fn), 1)
    )
    if denom == 0:
        return 0.0
    return (tp * tn - fp * fn) / denom


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", nargs="+", required=True)
    ap.add_argument("--results-dir", required=True,
                    help="Root of per-sample tool predictions "
                         "(e.g. results/test_real/full_run)")
    ap.add_argument("--gt-root", required=True,
                    help="Root of per-sample ground-truth dirs "
                         "(e.g. results/test_real/ground_truth)")
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--min-tier", type=int, default=2,
                    help="Minimum tier for viral positives (default 2)")
    ap.add_argument("--gt-variants", nargs="+",
                    default=["ground_truth.tsv", "ground_truth_no_e5.tsv"],
                    help="Per-sample GT filenames to compare")
    args = ap.parse_args()

    results_root = Path(args.results_dir)
    gt_root = Path(args.gt_root)
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    # tool -> variant -> [TP, FP, TN, FN] accumulator
    accum: dict[str, dict[str, list[int]]] = {}
    per_sample_rows: list[dict] = []

    for sample in args.samples:
        print(f"\n=== {sample} ===")
        # Load GT variants for this sample
        gts: dict[str, dict] = {}
        for variant in args.gt_variants:
            gt_path = gt_root / sample / variant
            if not gt_path.exists():
                print(f"  [WARN] missing {gt_path}; skipping this variant")
                continue
            gts[variant] = load_ground_truth(str(gt_path), min_tier=args.min_tier)
            n_viral = sum(1 for v in gts[variant].values() if v["label"] == "viral")
            n_neg = sum(1 for v in gts[variant].values() if v["label"] == "negative")
            print(f"  {variant:30s}  viral={n_viral:4d}  neg={n_neg:5d}")

        if not gts:
            continue

        # Union of GT contig IDs across variants (for tool parsing)
        any_variant = next(iter(gts.values()))
        gt_ids = set(any_variant.keys())
        for g in gts.values():
            gt_ids |= set(g.keys())

        # Parse tool predictions once for this sample
        # (parsers expect results_dir to be the per-sample subdir)
        sample_dir = results_root / sample
        tool_scores = discover_tool_outputs(
            sample_dir, sample, gt_ids=gt_ids, tools=None
        )
        if not tool_scores:
            print(f"  [WARN] no tool predictions found in {results_root}/{sample}")
            continue

        # Fallback for Jaeger: parser uses non-recursive glob but outputs live
        # at jaeger/{sample}_contigs/*_default_jaeger.tsv
        if "Jaeger" not in tool_scores or not tool_scores.get("Jaeger"):
            jaeger_scores = _parse_jaeger_recursive(sample_dir / "jaeger")
            if jaeger_scores:
                tool_scores["Jaeger"] = jaeger_scores
                print(f"  Jaeger (rglob fallback):       {len(jaeger_scores):5d} contigs parsed")

        for tool, scores in tool_scores.items():
            thr = TOOL_THRESHOLDS.get(tool, 0.5)
            pvalue_thr = TOOL_PVALUE_THRESHOLDS.get(tool)
            accum.setdefault(tool, {})
            for variant, gt in gts.items():
                tp, fp, tn, fn = confusion(scores, gt, thr, pvalue_thr)
                accum[tool].setdefault(variant, [0, 0, 0, 0])
                accum[tool][variant][0] += tp
                accum[tool][variant][1] += fp
                accum[tool][variant][2] += tn
                accum[tool][variant][3] += fn
                per_sample_rows.append({
                    "sample": sample, "tool": tool, "variant": variant,
                    "threshold": thr,
                    "TP": tp, "FP": fp, "TN": tn, "FN": fn,
                    "MCC": round(mcc(tp, fp, tn, fn), 4),
                })

    # Save per-sample long-format table
    per_sample_df = pd.DataFrame(per_sample_rows)
    per_sample_df.to_csv(outdir / "per_sample_confusion.tsv", sep="\t", index=False)

    # Pool confusion matrices and compute pooled MCC / rank per variant
    pooled_rows = []
    for tool, by_variant in accum.items():
        for variant, (tp, fp, tn, fn) in by_variant.items():
            pooled_rows.append({
                "tool": tool, "variant": variant,
                "TP": tp, "FP": fp, "TN": tn, "FN": fn,
                "precision": round(tp / max(tp + fp, 1), 4),
                "recall": round(tp / max(tp + fn, 1), 4),
                "MCC": round(mcc(tp, fp, tn, fn), 4),
            })
    pooled_df = pd.DataFrame(pooled_rows)
    pooled_df["rank"] = pooled_df.groupby("variant")["MCC"] \
        .rank(method="min", ascending=False).astype(int)
    pooled_df = pooled_df.sort_values(["variant", "rank"])
    pooled_df.to_csv(outdir / "pooled_by_variant.tsv", sep="\t", index=False)

    # Wide summary: one row per tool, MCC + rank under each variant
    wide = pooled_df.pivot_table(
        index="tool", columns="variant", values=["MCC", "rank"]
    )
    wide.columns = [f"{metric}_{variant.replace('ground_truth', 'GT').replace('.tsv','')}"
                    for metric, variant in wide.columns]
    variants = args.gt_variants
    if len(variants) == 2:
        col_primary = f"MCC_{variants[0].replace('ground_truth','GT').replace('.tsv','')}"
        col_no_e5 = f"MCC_{variants[1].replace('ground_truth','GT').replace('.tsv','')}"
        rank_primary = f"rank_{variants[0].replace('ground_truth','GT').replace('.tsv','')}"
        rank_no_e5 = f"rank_{variants[1].replace('ground_truth','GT').replace('.tsv','')}"
        wide["delta_MCC"] = (wide[col_no_e5] - wide[col_primary]).round(4)
        wide["delta_rank"] = wide[rank_no_e5] - wide[rank_primary]
        wide = wide.sort_values(col_primary, ascending=False)
    wide.to_csv(outdir / "sensitivity_s2_no_e5_summary.tsv", sep="\t")

    # Print compact summary to stdout
    print("\n" + "=" * 70)
    print("Pooled summary (MCC and rank) per tool per GT variant:")
    print("=" * 70)
    print(wide.to_string())
    print("\nOutputs written to", outdir)


if __name__ == "__main__":
    main()
