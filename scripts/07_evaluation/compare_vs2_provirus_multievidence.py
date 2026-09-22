#!/usr/bin/env python3
"""
compare_vs2_provirus_multievidence.py

Head-to-head comparison of VirSorter2 on a per-sample multi-evidence
benchmark contig set under two settings:

  - "off"  : original published runs with --provirus-off (provirus disabled)
  - "on"   : sensitivity re-run at developer default (provirus enabled)

Ground truth comes from the multi-evidence tier table (NOT FASTA labels);
the H2 hypothesis lives here, so the script reports recall on Tier 1
prophage contigs explicitly, alongside overall MCC.

Usage:
    python compare_vs2_provirus_multievidence.py \
        --sample UC093_V3 \
        --gt    results/test_real/ground_truth/UC093_V3/ground_truth.tsv \
        --off-tsv results/test_real/full_run/UC093_V3/virsorter2/final-viral-score.tsv \
        --on-tsv  results/test_real/full_run/UC093_V3/virsorter2_provirus_on/final-viral-score.tsv

Prints a Markdown-style summary to stdout.
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path


def parse_ground_truth(tsv: Path):
    """Return dict[contig_id] -> dict(tier=int, category=str)."""
    gt = {}
    with open(tsv) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            i_id   = header.index("contig_id")
            i_tier = header.index("tier")
            i_cat  = header.index("category")
        except ValueError as exc:
            raise SystemExit(f"ground_truth.tsv missing column: {exc}")
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            try:
                tier = int(cols[i_tier])
            except ValueError:
                continue
            gt[cols[i_id]] = {"tier": tier, "category": cols[i_cat]}
    return gt


def parse_vs2_scores(tsv: Path):
    """Return {contig_id: max_score}. Strips VS2's '||full' / '||N_partial' suffix."""
    scores = {}
    with open(tsv) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            i_name = header.index("seqname")
            i_max  = header.index("max_score")
        except ValueError as exc:
            raise SystemExit(f"VS2 TSV missing column: {exc}")
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            base = cols[i_name].split("||")[0]
            try:
                s = float(cols[i_max])
            except ValueError:
                continue
            if s > scores.get(base, -1.0):
                scores[base] = s
    return scores


def metrics(tp, fp, tn, fn):
    prec = tp / (tp + fp) if (tp + fp) else 0.0
    rec  = tp / (tp + fn) if (tp + fn) else 0.0
    f1   = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
    denom = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc = ((tp * tn) - (fp * fn)) / denom if denom else 0.0
    return prec, rec, f1, mcc


def confusion_against_tiers(gt, scores, threshold=0.5, positive_tiers=(1, 2),
                            positive_categories=("prophage", "free_phage", "eukaryotic_virus")):
    """
    TP = ground-truth Tier 1/2 viral contig called viral by VS2
    FN = ground-truth Tier 1/2 viral contig NOT called by VS2
    FP = ground-truth Tier 0 bacterial called viral by VS2
    TN = ground-truth Tier 0 bacterial NOT called by VS2

    Tier 3 (single-evidence) and Tier -1 (ambiguous) are excluded — that's
    how the multi-evidence benchmark in the manuscript treats them too.
    """
    tp = fp = tn = fn = 0
    for contig, info in gt.items():
        called = scores.get(contig, 0.0) >= threshold
        is_pos = info["tier"] in positive_tiers and info["category"] in positive_categories
        is_neg = info["tier"] == 0 and info["category"] == "bacterial"
        if is_pos and called:
            tp += 1
        elif is_pos and not called:
            fn += 1
        elif is_neg and called:
            fp += 1
        elif is_neg and not called:
            tn += 1
        # else: excluded (Tier 3, Tier -1)
    return tp, fp, tn, fn


def tier1_prophage_recall(gt, scores, threshold=0.5):
    """The exact H2 metric: VS2 recall on Tier 1 prophage contigs."""
    targets = [c for c, info in gt.items()
               if info["tier"] == 1 and info["category"] == "prophage"]
    hit = sum(1 for c in targets if scores.get(c, 0.0) >= threshold)
    return hit, len(targets)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample",  required=True)
    ap.add_argument("--gt",      required=True, type=Path)
    ap.add_argument("--off-tsv", required=True, type=Path)
    ap.add_argument("--on-tsv",  required=True, type=Path)
    ap.add_argument("--threshold", type=float, default=0.5)
    args = ap.parse_args()

    gt  = parse_ground_truth(args.gt)
    off = parse_vs2_scores(args.off_tsv)
    on  = parse_vs2_scores(args.on_tsv)

    n_total = len(gt)
    n_tier1_pro = sum(1 for v in gt.values() if v["tier"] == 1 and v["category"] == "prophage")
    n_tier1_fp  = sum(1 for v in gt.values() if v["tier"] == 1 and v["category"] == "free_phage")
    n_tier1_euk = sum(1 for v in gt.values() if v["tier"] == 1 and v["category"] == "eukaryotic_virus")
    n_tier0_bact= sum(1 for v in gt.values() if v["tier"] == 0 and v["category"] == "bacterial")

    print(f"# VirSorter2 provirus-mode sensitivity — sample {args.sample}\n")
    print(f"- Ground-truth contigs        : {n_total}")
    print(f"- Tier 1 prophage (H2 target) : {n_tier1_pro}")
    print(f"- Tier 1 free_phage           : {n_tier1_fp}")
    print(f"- Tier 1 eukaryotic_virus     : {n_tier1_euk}")
    print(f"- Tier 0 bacterial (TN pool)  : {n_tier0_bact}")
    print(f"- Threshold (max_score ≥)     : {args.threshold}\n")

    print("## H2 metric — Tier 1 prophage recall (the headline number)\n")
    h_off, _ = tier1_prophage_recall(gt, off, args.threshold)
    h_on,  n = tier1_prophage_recall(gt, on,  args.threshold)
    delta = h_on - h_off
    print(f"| Setting | Recall on Tier 1 prophage |")
    print(f"|---|---:|")
    print(f"| --provirus-off (published)        | {h_off} / {n} = {h_off/max(n,1)*100:.1f}% |")
    print(f"| --provirus default (sensitivity)  | {h_on} / {n} = {h_on/max(n,1)*100:.1f}% |")
    print(f"| **Δ (sensitivity − published)**   | **{delta:+d}** contigs |\n")

    print("## Overall classification (Tier 1+2 viral vs Tier 0 bacterial)\n")
    print(f"| Setting | TP | FP | TN | FN | Precision | Recall | F1 | MCC |")
    print(f"|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    for name, sc in [("off (published)", off), ("on (sensitivity)", on)]:
        tp, fp, tn, fn = confusion_against_tiers(gt, sc, args.threshold)
        p, r, f1, mcc = metrics(tp, fp, tn, fn)
        print(f"| {name} | {tp} | {fp} | {tn} | {fn} | {p:.3f} | {r:.3f} | {f1:.3f} | {mcc:.3f} |")
    print()

    # Per-contig flow: who moved between the two settings
    off_called = {c for c, s in off.items() if s >= args.threshold}
    on_called  = {c for c, s in on.items()  if s >= args.threshold}
    only_on  = on_called  - off_called
    only_off = off_called - on_called
    print("## Per-contig flow between the two settings\n")
    print(f"- Called viral in BOTH      : {len(off_called & on_called)}")
    print(f"- Called viral ONLY in OFF  : {len(only_off)}")
    print(f"- Called viral ONLY in ON   : {len(only_on)}")
    print()

    def stratify(s, label):
        if not s:
            return
        print(f"### Stratification of '{label}' by ground truth\n")
        bins = {}
        for c in s:
            info = gt.get(c)
            if info is None:
                key = ("not-in-GT", "—")
            else:
                key = (f"tier{info['tier']}", info["category"])
            bins[key] = bins.get(key, 0) + 1
        for (tier, cat), n_ in sorted(bins.items(), key=lambda kv: -kv[1]):
            print(f"- {tier} / {cat}: {n_}")
        print()

    stratify(only_on,  f"newly called viral by enabling provirus ({len(only_on)})")
    stratify(only_off, f"called viral by --provirus-off but NOT by default ({len(only_off)})")


if __name__ == "__main__":
    main()
