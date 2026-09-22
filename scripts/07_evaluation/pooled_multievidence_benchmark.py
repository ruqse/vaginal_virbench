#!/usr/bin/env python3
"""Pooled multi-evidence benchmark across all 13 matched-virome samples.

Replaces the original single-sample (UC028_V2) "secondary benchmark" with a
genuine 13-sample pooled evaluation, so the manuscript's "pooled across all 13
participants" framing becomes true.

Design (avoids the lockstep contig-id namespacing hazard):
  * Metrics are pooled by SUMMING per-sample confusion counts and by appending
    per-sample (label, score) to pooled vectors for AUPRC. There is no global
    contig_id dict, so cross-sample id collisions (144 bare NODE_* ids collide)
    cannot silently drop or mismatch anything.
  * Contig ids are namespaced "<sample>|<contig_id>" ONLY when writing the audit
    artifacts, never during matching.

Reuses parsers/metrics from evaluate_metagenome.py and the recursive Jaeger
parser + _predict from sensitivity_no_e5.py.

Gate 0 (length policy) is an explicit --min-length parameter recorded in REPORT.md.

Usage:
    module load SciPy-bundle/2024.05-gfbf-2024a
    python scripts/07_evaluation/pooled_multievidence_benchmark.py \
        --samples UC028_V2 UC055_V1 UC055_V2 UC062_V2 UC065_V2 UC074_V2 \
                  UC084_V2 UC093_V2 UC093_V3 UC096_V2 UC115_V2 UC139_V2 UC164_V2 \
        --results-dir results/test_real/full_run \
        --gt-root     results/test_real/ground_truth \
        --min-length  1500 \
        --output-dir  results/test_real/secondary_benchmark_pooled13
"""
from __future__ import annotations

import argparse
import csv
import gzip
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from evaluate_metagenome import (  # noqa: E402
    TOOL_PVALUE_THRESHOLDS,
    TOOL_THRESHOLDS,
    _auprc,
    discover_tool_outputs,
    load_ground_truth,
)
from sensitivity_no_e5 import _parse_jaeger_recursive, _predict  # noqa: E402
from deterministic_gzip import DeterministicGzipWriter  # noqa: E402

ALL_TOOLS = list(TOOL_THRESHOLDS.keys())  # all 14, even if a tool parsed nothing


def _mcc(tp: int, fp: int, tn: int, fn: int) -> float:
    denom = math.sqrt(
        max(tp + fp, 1) * max(tp + fn, 1) * max(tn + fp, 1) * max(tn + fn, 1)
    )
    return (tp * tn - fp * fn) / denom if denom > 0 else 0.0


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", nargs="+", required=True)
    ap.add_argument("--results-dir", default="results/test_real/full_run")
    ap.add_argument("--gt-root", default="results/test_real/ground_truth")
    ap.add_argument("--gt-name", default="ground_truth_with_kraken2.tsv")
    ap.add_argument("--min-length", type=int, default=1500,
                    help="GATE 0: minimum contig length to include (default 1500 bp, "
                         "honouring Methods L131; use 0 for the original no-filter policy)")
    ap.add_argument("--min-tier", type=int, default=2)
    ap.add_argument("--output-dir",
                    default="results/test_real/secondary_benchmark_pooled13")
    args = ap.parse_args()

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    accum = {t: [0, 0, 0, 0] for t in ALL_TOOLS}      # tool -> [TP,FP,TN,FN]
    yt = {t: [] for t in ALL_TOOLS}                   # pooled AUPRC labels
    ys = {t: [] for t in ALL_TOOLS}                   # pooled AUPRC scores
    matched = {t: 0 for t in ALL_TOOLS}               # contigs with an emitted call
    gt_rows = []                                       # benchmark-ready GT artifact
    n_total = 0
    n_viral = n_neg = n_t1 = n_t2 = 0
    per_sample = []                                    # per-sample audit summary

    pred_path = outdir / "pooled_tool_predictions.tsv.gz"
    pred_fh = DeterministicGzipWriter(pred_path, newline="")  # byte-reproducible .gz
    pred_w = csv.writer(pred_fh, delimiter="\t")
    pred_w.writerow(["sample", "contig_id", "orig_contig_id", "tool", "label",
                     "score", "pvalue", "threshold", "y_pred", "emitted_call"])

    for sample in args.samples:
        gt_path = Path(args.gt_root) / sample / args.gt_name
        if not gt_path.exists():
            print(f"[FATAL] missing GT: {gt_path}", file=sys.stderr)
            sys.exit(1)
        gt = load_ground_truth(str(gt_path), min_tier=args.min_tier)
        # GATE 0 length filter
        gt = {c: v for c, v in gt.items() if v["length"] >= args.min_length}
        gt_ids = set(gt.keys())
        sv = sum(1 for v in gt.values() if v["label"] == "viral")
        sn = sum(1 for v in gt.values() if v["label"] == "negative")
        print(f"\n=== {sample}: {len(gt)} contigs (viral={sv}, neg={sn}) "
              f"@ >={args.min_length}bp ===")

        sample_dir = Path(args.results_dir) / sample
        tool_scores = discover_tool_outputs(sample_dir, sample, gt_ids=gt_ids)
        if not tool_scores.get("Jaeger"):
            jg = _parse_jaeger_recursive(sample_dir / "jaeger")
            if jg:
                tool_scores["Jaeger"] = jg
                print(f"  Jaeger (rglob fallback):       {len(jg):6d} contigs parsed")

        # GT artifact rows + global counts
        for c, info in gt.items():
            gt_rows.append({
                "sample": sample, "contig_id": f"{sample}|{c}", "orig_contig_id": c,
                "tier": info["tier"], "category": info["category"],
                "label": info["label"], "length": info["length"],
                "evidence_lines": info.get("evidence_lines", "-"),
            })
        n_total += len(gt)
        n_viral += sv
        n_neg += sn
        n_t1 += sum(1 for v in gt.values() if v["tier"] == 1)
        n_t2 += sum(1 for v in gt.values() if v["tier"] == 2)

        s_match = {}
        for tool in ALL_TOOLS:
            scores = tool_scores.get(tool, {})
            thr = TOOL_THRESHOLDS[tool]
            pvthr = TOOL_PVALUE_THRESHOLDS.get(tool)
            tp = fp = tn = fn = 0
            m = 0
            for c, info in gt.items():
                emitted = c in scores
                sc = scores.get(c, 0.0)
                s_comp = float(sc[0]) if isinstance(sc, tuple) else float(sc)
                pv = sc[1] if isinstance(sc, tuple) else ""
                pred = _predict(sc, thr, pvthr)
                yt[tool].append(1 if info["label"] == "viral" else 0)
                ys[tool].append(s_comp)
                if emitted:
                    m += 1
                if info["label"] == "viral":
                    tp += pred
                    fn += not pred
                else:
                    fp += pred
                    tn += not pred
                pred_w.writerow([sample, f"{sample}|{c}", c, tool, info["label"],
                                 f"{s_comp:.6g}", pv, thr, int(bool(pred)), int(emitted)])
            accum[tool][0] += tp
            accum[tool][1] += fp
            accum[tool][2] += tn
            accum[tool][3] += fn
            matched[tool] += m
            s_match[tool] = m
        per_sample.append({"sample": sample, "n": len(gt), "viral": sv, "neg": sn,
                           **{f"matched_{t}": s_match[t] for t in ALL_TOOLS}})

    pred_fh.close()

    # --- Pooled per-tool metrics ---
    overall = []
    for tool in ALL_TOOLS:
        tp, fp, tn, fn = accum[tool]
        assert tp + fp + tn + fn == n_total, (
            f"{tool}: confusion sum {tp+fp+tn+fn} != n_total {n_total}")
        prec = tp / (tp + fp) if (tp + fp) else 0.0
        rec = tp / (tp + fn) if (tp + fn) else 0.0
        f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
        yt_a = np.array(yt[tool])
        ys_a = np.array(ys[tool])
        auprc = round(_auprc(yt_a, ys_a), 4) if len(np.unique(ys_a)) > 1 else "NA"
        overall.append({
            "tool": tool, "TP": tp, "FP": fp, "TN": tn, "FN": fn,
            "precision": round(prec, 4), "recall": round(rec, 4),
            "F1": round(f1, 4), "MCC": round(_mcc(tp, fp, tn, fn), 4),
            "AUPRC": auprc, "n_predicted": tp + fp, "n_total": n_total,
            "n_emitted": matched[tool],
        })
    overall.sort(key=lambda r: r["MCC"], reverse=True)

    # --- Write outputs ---
    def write_tsv(rows, path, cols):
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
            w.writeheader()
            w.writerows(rows)

    write_tsv(overall, outdir / "overall_metrics.tsv",
              ["tool", "TP", "FP", "TN", "FN", "precision", "recall", "F1",
               "MCC", "AUPRC", "n_predicted", "n_total", "n_emitted"])
    write_tsv(gt_rows, outdir / "pooled_benchmark_ready_ground_truth.tsv",
              ["sample", "contig_id", "orig_contig_id", "tier", "category",
               "label", "length", "evidence_lines"])
    write_tsv(per_sample, outdir / "per_sample_summary.tsv",
              ["sample", "n", "viral", "neg"] + [f"matched_{t}" for t in ALL_TOOLS])

    ratio = n_neg / max(n_viral, 1)
    with open(outdir / "REPORT.md", "w") as fh:
        fh.write("# Pooled multi-evidence benchmark (13 samples)\n\n")
        fh.write("**GATE 0 — contig length policy**: "
                 f"`--min-length {args.min_length}` "
                 f"({'>=%d bp (Methods L131)' % args.min_length if args.min_length else 'no filter'}).\n\n")
        fh.write(f"- Samples: {len(args.samples)} ({', '.join(args.samples)})\n")
        fh.write(f"- min_tier (viral positives): tier 1..{args.min_tier}\n")
        fh.write(f"- **Total benchmark-ready contigs**: {n_total}\n")
        fh.write(f"- **Viral positives (Tier 1+2)**: {n_viral} "
                 f"(Tier 1 = {n_t1}, Tier 2 = {n_t2})\n")
        fh.write(f"- **True negatives (Tier 0)**: {n_neg}\n")
        fh.write(f"- **Class ratio**: 1:{ratio:.1f} (viral:negative)\n\n")
        fh.write("## Overall metrics (ordered by MCC)\n\n")
        fh.write("| Tool | TP | FP | TN | FN | Prec | Rec | F1 | MCC | AUPRC | n_emitted |\n")
        fh.write("|---|---|---|---|---|---|---|---|---|---|---|\n")
        for r in overall:
            fh.write(f"| {r['tool']} | {r['TP']} | {r['FP']} | {r['TN']} | {r['FN']} "
                     f"| {r['precision']} | {r['recall']} | {r['F1']} | {r['MCC']} "
                     f"| {r['AUPRC']} | {r['n_emitted']} |\n")

    print(f"\n{'=' * 70}")
    print(f"Pooled benchmark complete: n={n_total} (viral={n_viral} "
          f"[T1={n_t1},T2={n_t2}], neg={n_neg}, 1:{ratio:.1f})")
    print(f"Top 4 by MCC: " + ", ".join(
        f"{r['tool']} {r['MCC']}" for r in overall[:4]))
    print(f"Outputs in: {outdir}")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
