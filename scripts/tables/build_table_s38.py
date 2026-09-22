#!/usr/bin/env python3
"""Build Table S38 -- per-tool false-positive rate by negative-grounding category.

Why this table exists
---------------------
Figure S14 lists "Kraken2 human" first in its negative-grounding box, and Figure
S15B shows it as the largest Kraken2 category among the multi-evidence negatives
(2,833 of 11,577, 24%; three times the 964 bacterial/archaeal). Neither figure,
nor any numbered table, reported what the tools actually DO on those contigs.
The manuscript separately states that human assignments do not qualify as Track C
negatives, so a reader meets human DNA as a first-class negative class in one
benchmark and as an excluded class in the other, with no per-tool numbers for
either. S38 supplies them.

The result is not a restatement of overall specificity. Tools separate far more
sharply on human contigs (1% to 96%) than the pooled negative FPR suggests, and
two tools inverts against their own bacterial performance: TransGINmer is much
worse on human than on bacterial contigs, while Jaeger is much better.

Created as Table S39; renumbered to S38 on 2026-09-15 when the former Table S38
(cross-study design audit) was removed from the supplement.

Sources, both already in the repository; no tool is re-run:
  results/test_real/ground_truth/<sample>/ground_truth_with_kraken2.tsv
      Per-contig tier and `notes`. The grounding branch order is copied verbatim
      from scripts/07_evaluation/viralm_investigation/06_plot_gt_construction.py
      (human -> kraken2: -> checkv) so this table reproduces Figure S15B exactly.
  results/test_real/secondary_benchmark_pooled13/pooled_tool_predictions.tsv.gz
      Complete 13,030 x 14 contig-by-tool matrix carrying `y_pred` already
      evaluated at the documented operating thresholds (Table S21). Using the
      stored y_pred rather than re-thresholding keeps S38 consistent with
      Table S12 by construction.

The Track C columns repeat the measurement on an independently labelled negative
set: the human-assigned contigs that fall in the Track C background enrichment
band. They are a replication check, not part of the multi-evidence benchmark, and
they are few (n = 80). Track C calls come from the Track C evaluator's own
parsers so the operating rules match that track.

    python3 scripts/tables/build_table_s38.py
"""
import csv
import gzip
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]         # repo root (scripts/tables/..)
HERE = ROOT / "results" / "tables"                 # published tables
sys.path.insert(0, str(ROOT / "scripts/03_assembly"))

GTROOT = ROOT / "results/test_real/ground_truth"
PRED = ROOT / "results/test_real/secondary_benchmark_pooled13/pooled_tool_predictions.tsv.gz"
TRACKC_GT = ROOT / "results/test_real/coassembly/rca_ground_truth.tsv"
TRACKC_RESULTS = ROOT / "results/test_real/full_run/RCA_master"
TAX = HERE / "tool_method_taxonomy_verified.tsv"
OUT = HERE / "table_s38_negative_grounding_specificity.tsv"

MIN_LEN = 1500
SAMPLES = [
    "UC028_V2", "UC055_V1", "UC055_V2", "UC062_V2", "UC065_V2",
    "UC074_V2", "UC084_V2", "UC093_V2", "UC093_V3", "UC096_V2",
    "UC115_V2", "UC139_V2", "UC164_V2",
]

# Display labels for the three grounding categories, in Figure S15B order.
CATS = [("human", "human"), ("bact_arch", "bacterial/archaeal"), ("checkv_host", "CheckV host genes")]


def negative_grounding() -> dict:
    """Return {pooled_contig_id: category} for the Tier 0 negatives.

    Branch order is copied from 06_plot_gt_construction.py so the category
    totals reproduce Figure S15B (2,833 / 964 / 7,780).
    """
    grounding = {}
    for sample in SAMPLES:
        path = GTROOT / sample / "ground_truth_with_kraken2.tsv"
        if not path.exists():
            raise FileNotFoundError(f"missing ground-truth source table: {path}")
        with path.open(newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if int(row["length"]) < MIN_LEN or int(row["tier"]) != 0:
                    continue
                notes = row.get("notes", "") or ""
                if "Homo sapiens" in notes or "9606" in notes:
                    cat = "human"
                elif notes.startswith("kraken2:"):
                    cat = "bact_arch"
                elif notes.startswith("checkv"):
                    cat = "checkv_host"
                else:
                    continue
                # The prediction matrix joins on "<sample>|<contig>".
                grounding[f"{sample}|{row['contig_id']}"] = cat
    return grounding


def multievidence_fp(grounding: dict):
    """Per-tool false positives on each grounding category, from stored y_pred."""
    fp = defaultdict(Counter)
    n = defaultdict(Counter)
    emitted = defaultdict(int)
    with gzip.open(PRED, "rt", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["label"] != "negative":
                continue
            cat = grounding.get(row["contig_id"])
            if cat is None:
                continue
            tool = row["tool"]
            n[tool][cat] += 1
            if row["y_pred"] == "1":
                fp[tool][cat] += 1
            if row["emitted_call"] == "1":
                emitted[tool] += 1
    return fp, n, emitted


def trackc_human_fp():
    """Replication on the Track C human contigs in the background band."""
    import evaluate_rca_benchmark as E

    human = set()
    with TRACKC_GT.open(newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if ("Homo sapiens" in row["kraken2_taxon"]
                    and row["enrichment_class"] == "background"
                    and float(row["viral_genes"] or 0) == 0):
                human.add(row["contig_id"])

    # The Track C tool outputs are not redistributed; without them this table
    # cannot be rebuilt, so stop instead of silently writing zero counts.
    if not TRACKC_RESULTS.is_dir():
        raise FileNotFoundError(f"Track C tool outputs not found: {TRACKC_RESULTS} "
                                "(Table S38 is shipped precomputed)")
    scores = E.discover_tool_outputs(TRACKC_RESULTS, "RCA_master", gt_ids=human)
    if not scores:
        raise FileNotFoundError(f"no Track C tool outputs found under {TRACKC_RESULTS}")
    out = {}
    for tool, sc in scores.items():
        thr = E.TOOL_THRESHOLDS.get(tool, 0.5)
        pthr = E.TOOL_PVALUE_THRESHOLDS.get(tool)
        hits = 0
        for cid in human:
            val = sc.get(cid)
            if val is None:
                continue
            if isinstance(val, tuple):
                score, pv = val[0], val[1]
                if score >= thr and (pthr is None or pv < pthr):
                    hits += 1
            elif val >= thr:
                hits += 1
        out[tool] = hits
    return out, len(human)


def approaches() -> dict:
    tax = {}
    with TAX.open(newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            tax[row["tool"]] = row["method_approach"]
    return tax


def rate(num, den):
    return f"{100 * num / den:.1f}" if den else "NA"


def main():
    grounding = negative_grounding()
    totals = Counter(grounding.values())
    print("negative grounding (must match Fig. S15B):", dict(totals))

    fp, n, emitted = multievidence_fp(grounding)
    tc_fp, tc_n = trackc_human_fp()
    tax = approaches()

    rows = []
    for tool in sorted(n, key=lambda t: -(fp[t]["human"] / max(n[t]["human"], 1))):
        h_fpr = 100 * fp[tool]["human"] / n[tool]["human"]
        b_fpr = 100 * fp[tool]["bact_arch"] / n[tool]["bact_arch"]
        rows.append({
            "tool": tool,
            "method_approach": tax.get(tool, "NA"),
            "contigs_reported_by_tool": emitted[tool],
            "human_n": n[tool]["human"],
            "human_fp": fp[tool]["human"],
            "human_fpr_pct": rate(fp[tool]["human"], n[tool]["human"]),
            "bact_arch_n": n[tool]["bact_arch"],
            "bact_arch_fp": fp[tool]["bact_arch"],
            "bact_arch_fpr_pct": rate(fp[tool]["bact_arch"], n[tool]["bact_arch"]),
            "checkv_host_n": n[tool]["checkv_host"],
            "checkv_host_fp": fp[tool]["checkv_host"],
            "checkv_host_fpr_pct": rate(fp[tool]["checkv_host"], n[tool]["checkv_host"]),
            "human_minus_bact_arch_pp": f"{h_fpr - b_fpr:+.1f}",
            "trackc_human_n": tc_n,
            "trackc_human_fp": tc_fp.get(tool, 0),
            "trackc_human_fpr_pct": rate(tc_fp.get(tool, 0), tc_n),
        })

    with OUT.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {OUT.relative_to(ROOT)} ({len(rows)} tools)")


if __name__ == "__main__":
    main()
