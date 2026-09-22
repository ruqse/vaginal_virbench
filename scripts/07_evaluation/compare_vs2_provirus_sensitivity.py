#!/usr/bin/env python3
"""
compare_vs2_provirus_sensitivity.py

Head-to-head comparison of VirSorter2 on the same input under two settings:
  - "off"  : original published runs with --provirus-off (provirus disabled)
  - "on"   : sensitivity re-run at developer default (provirus enabled)

Reports per-fragment agreement plus standard classification metrics (TP/FP/TN/
FN, precision, recall, F1, MCC) against the ground-truth labels encoded in the
fragment FASTA headers (label=viral / label=non-viral).

Usage:
    python compare_vs2_provirus_sensitivity.py \
        --fasta data/spike_in/fragments/L1500_fragments.fasta \
        --off-tsv results/spike_in_benchmark/L1500_fragments/virsorter2/final-viral-score.tsv \
        --on-tsv  results/spike_in_benchmark/L1500_fragments/virsorter2_provirus_on/final-viral-score.tsv

Prints a Markdown-style summary to stdout.
"""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path


def parse_fasta_labels(fasta: Path) -> dict[str, str]:
    """Return {fragment_id: 'viral'|'non-viral'} from FASTA headers."""
    labels: dict[str, str] = {}
    label_re = re.compile(r"label=(\S+)")
    with open(fasta) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            head = line[1:].rstrip("\n")
            frag_id = head.split()[0]
            m = label_re.search(head)
            if not m:
                raise SystemExit(f"No label= tag on header: {head!r}")
            labels[frag_id] = m.group(1)
    return labels


def parse_vs2_scores(tsv: Path) -> dict[str, float]:
    """Return {fragment_id: max_score} for all fragments VS2 called viral.

    The 'seqname' field carries a '||full' or '||{n}_partial' suffix from
    VS2. We strip everything after '||' to recover the base fragment ID and
    take the max score across any rows that map to the same fragment.
    """
    scores: dict[str, float] = {}
    with open(tsv) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            i_name = header.index("seqname")
            i_max  = header.index("max_score")
        except ValueError as exc:
            raise SystemExit(f"VS2 TSV missing expected column: {exc}")
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


def confusion(labels: dict[str, str], calls: dict[str, float], threshold: float = 0.5):
    tp = fp = tn = fn = 0
    for frag, gt in labels.items():
        called = calls.get(frag, 0.0) >= threshold
        if gt == "viral" and called:
            tp += 1
        elif gt == "viral" and not called:
            fn += 1
        elif gt != "viral" and called:
            fp += 1
        else:
            tn += 1
    return tp, fp, tn, fn


def metrics(tp: int, fp: int, tn: int, fn: int):
    prec = tp / (tp + fp) if (tp + fp) else 0.0
    rec  = tp / (tp + fn) if (tp + fn) else 0.0
    f1   = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
    denom = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc = ((tp * tn) - (fp * fn)) / denom if denom else 0.0
    return prec, rec, f1, mcc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True, type=Path)
    ap.add_argument("--off-tsv", required=True, type=Path,
                    help="VirSorter2 final-viral-score.tsv from --provirus-off run")
    ap.add_argument("--on-tsv", required=True, type=Path,
                    help="VirSorter2 final-viral-score.tsv from default (provirus on) run")
    ap.add_argument("--threshold", type=float, default=0.5,
                    help="max_score cutoff (default 0.5, matches manuscript)")
    args = ap.parse_args()

    labels = parse_fasta_labels(args.fasta)
    off    = parse_vs2_scores(args.off_tsv)
    on     = parse_vs2_scores(args.on_tsv)

    n_viral     = sum(1 for v in labels.values() if v == "viral")
    n_non_viral = sum(1 for v in labels.values() if v != "viral")

    off_called = {f for f, s in off.items() if s >= args.threshold and f in labels}
    on_called  = {f for f, s in on.items()  if s >= args.threshold and f in labels}

    only_off = off_called - on_called
    only_on  = on_called  - off_called
    both     = off_called & on_called

    print("# VirSorter2 provirus-mode sensitivity comparison\n")
    print(f"- Input FASTA       : {args.fasta}")
    print(f"- Off-mode TSV      : {args.off_tsv}")
    print(f"- On-mode  TSV      : {args.on_tsv}")
    print(f"- Threshold (max_score ≥) : {args.threshold}\n")

    print(f"- Total fragments      : {len(labels):>6} ({n_viral} viral, {n_non_viral} non-viral)\n")

    print("## Calls per setting\n")
    print(f"| Setting | Viral calls | of which TP | of which FP |")
    print(f"|---|---:|---:|---:|")
    for name, called in [("--provirus-off (published)", off_called),
                         ("--provirus default (sensitivity)", on_called)]:
        tp = sum(1 for f in called if labels[f] == "viral")
        fp = len(called) - tp
        print(f"| {name} | {len(called)} | {tp} | {fp} |")
    print()

    print("## Per-fragment agreement\n")
    print(f"- Called viral in BOTH      : {len(both):>5}")
    print(f"- Called viral ONLY in OFF  : {len(only_off):>5}")
    print(f"- Called viral ONLY in ON   : {len(only_on):>5}")
    print()

    print("## Classification metrics vs ground truth\n")
    print(f"| Setting | TP | FP | TN | FN | Precision | Recall | F1 | MCC |")
    print(f"|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    for name, calls in [("off (published)", off), ("on (sensitivity)", on)]:
        tp, fp, tn, fn = confusion(labels, calls, args.threshold)
        p, r, f1, mcc = metrics(tp, fp, tn, fn)
        print(f"| {name} | {tp} | {fp} | {tn} | {fn} | {p:.3f} | {r:.3f} | {f1:.3f} | {mcc:.3f} |")
    print()

    # Stratify the only-ON gains by ground truth — likely most are viral
    # fragments rescued by the provirus boundary classifier.
    if only_on:
        only_on_viral = sum(1 for f in only_on if labels[f] == "viral")
        only_on_bact  = len(only_on) - only_on_viral
        print(f"## Calls gained by enabling provirus detection\n")
        print(f"- {len(only_on)} fragments newly called viral with --provirus default:")
        print(f"  - {only_on_viral} truly viral (correct gains)")
        print(f"  - {only_on_bact} truly non-viral (false-positive gains)")
        print()
    if only_off:
        only_off_viral = sum(1 for f in only_off if labels[f] == "viral")
        only_off_bact  = len(only_off) - only_off_viral
        print(f"## Calls dropped by enabling provirus detection\n")
        print(f"- {len(only_off)} fragments called viral by --provirus-off but NOT default:")
        print(f"  - {only_off_viral} truly viral (regressions)")
        print(f"  - {only_off_bact} truly non-viral (FP that disappeared)")


if __name__ == "__main__":
    main()
