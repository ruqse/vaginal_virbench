#!/usr/bin/env python3
"""
06_audit_e1b_parser.py — narrow audit of the legacy Evidence-1b parser.

The legacy E1b parser (run_blastn_imgvr.sh) sums HSP lengths across DIFFERENT
subjects and never merges overlapping query intervals, so its aligned fraction
can be inflated -> a contig that should fail AF>=75% may pass. This audit
re-parses the EXISTING E1b BLAST output with the corrected interval-merged
parser and asks whether any GROUND-TRUTH TIER changes.

CRITICAL — model the real tier logic, do not treat E1b as a standalone line.
build_ground_truth.py collapses homology into a single line E1 = E1a OR E1b
(see assign_tier_and_category). So a changed E1b call only matters when it flips
has_e1 = has_e1a OR has_e1b -- it cannot for any contig E1a already supports.
We therefore report n_e1b_changed AND n_tier_changed (Tier-1 flips separately).

READ-ONLY: this script never re-runs BLAST and never writes ground_truth.tsv.

Usage:
  python 06_audit_e1b_parser.py \
      --gt-root results/test_real/ground_truth \
      --out results/vmgc_crosscheck/e1b_parser_audit.tsv
"""

import argparse
import csv
import glob
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "04_ground_truth"))

from vmgc_coverage import best_subject_per_query  # noqa: E402
import build_ground_truth as bgt  # noqa: E402

MIN_ANI = 90.0   # matches run_blastn_imgvr.sh
MIN_AF = 75.0


def corrected_e1b_membership(raw_blast):
    """Return set of contig_ids passing E1b under the corrected parser."""
    members = {}
    if not os.path.isfile(raw_blast):
        return members
    best = best_subject_per_query(raw_blast)
    for cid, m in best.items():
        if m["weighted_identity"] >= MIN_ANI and m["q_af"] >= MIN_AF:
            members[cid] = {"metavr_hit": m["sseqid"], "ani": m["weighted_identity"]}
    return members


def tier_of(cid, e1a, e1b, e2, e3, e4, e5, hg, kr):
    return bgt.assign_tier_and_category(
        cid, e1a, e1b, e2, e3, e4, e5, hg, kraken2_taxonomy=kr)[0]


def load_frozen_tiers(sample_dir):
    """contig_id -> int tier, from the primary frozen ground_truth.tsv."""
    gt = os.path.join(sample_dir, "ground_truth.tsv")
    tiers = {}
    if not os.path.isfile(gt):
        return tiers
    with open(gt) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            try:
                tiers[row["contig_id"]] = int(row["tier"])
            except (ValueError, KeyError):
                continue
    return tiers


def detect_e5(sample_dir):
    """Reproduce the PRIMARY ground truth's E5 source.
    The primary ground_truth.tsv was built from read-level E5 where only that
    file exists; prefer evidence_5_rca_readlevel.tsv, else assembly."""
    rl = os.path.join(sample_dir, "evidence_5_rca_readlevel.tsv")
    asm = os.path.join(sample_dir, "evidence_5_rca.tsv")
    if os.path.isfile(rl):
        return bgt.load_evidence_5(rl, mode="read")
    return bgt.load_evidence_5(asm, mode="assembly")


def audit_sample(sample_dir):
    raw = os.path.join(sample_dir, "blastn_imgvr_raw.tsv")
    if not os.path.isfile(raw):
        return None
    sample = os.path.basename(sample_dir.rstrip("/"))

    e1a = bgt.load_evidence_1a(os.path.join(sample_dir, "evidence_1a_diamond.tsv"))
    e1b_orig = bgt.load_evidence_1b(os.path.join(sample_dir, "evidence_1b_imgvr.tsv"))
    e2 = bgt.load_evidence_2(os.path.join(sample_dir, "evidence_2_structural.tsv"))
    e3 = bgt.load_evidence_3(os.path.join(sample_dir, "evidence_3_prophage.tsv"))
    e4 = bgt.load_evidence_4(os.path.join(sample_dir, "evidence_4_crispr.tsv"))
    e5 = detect_e5(sample_dir)  # read-level E5 to match the primary ground truth
    hg = bgt.load_checkv_host_genes(os.path.join(sample_dir, "evidence_2_structural.tsv"))
    # Primary ground_truth.tsv used --kraken2-taxonomy when kraken2_output.tsv exists
    # (downstream_when_ready.sh). Load it so the baseline reproduces the frozen GT.
    kr_path = os.path.join(sample_dir, "kraken2_output.tsv")
    kr = bgt.load_kraken2_taxonomy(kr_path) if os.path.isfile(kr_path) else {}
    frozen = load_frozen_tiers(sample_dir)

    e1b_corr = corrected_e1b_membership(raw)
    orig_set, corr_set = set(e1b_orig), set(e1b_corr)
    gained = corr_set - orig_set   # corrected adds E1b support
    lost = orig_set - corr_set     # corrected removes E1b support

    # Faithfulness: recomputed-original tier must reproduce the frozen tier.
    # The pure-E1b DELTA below is machinery-independent (orig vs corr, same code),
    # but we report concordance so the baseline is demonstrably faithful.
    universe = set(frozen) | orig_set | corr_set | set(e1a) | set(e2) | set(e3) | set(e4) | set(e5)
    n_faithful = n_checked = 0
    tier_changed = tier1_changed = 0
    detail = []
    for cid in universe:
        t_orig = tier_of(cid, e1a, e1b_orig, e2, e3, e4, e5, hg, kr)
        t_corr = tier_of(cid, e1a, e1b_corr, e2, e3, e4, e5, hg, kr)
        if cid in frozen:
            n_checked += 1
            n_faithful += int(t_orig == frozen[cid])
        if t_orig != t_corr:
            tier_changed += 1
            if (t_orig == 1) != (t_corr == 1):
                tier1_changed += 1
            detail.append({
                "sample": sample, "contig_id": cid,
                "e1b_change": "gained" if cid in gained else "lost",
                "has_e1a": cid in e1a,
                "tier_frozen": frozen.get(cid, "NA"),
                "tier_recomp_orig": t_orig, "tier_corrected": t_corr,
            })

    return {
        "sample": sample,
        "n_e1b_orig": len(orig_set),
        "n_e1b_corrected": len(corr_set),
        "n_e1b_gained": len(gained),
        "n_e1b_lost": len(lost),
        "n_e1b_changed": len(gained | lost),
        "frozen_concordance_pct": round(100.0 * n_faithful / n_checked, 3) if n_checked else 0.0,
        "n_tier_changed": tier_changed,
        "n_tier1_changed": tier1_changed,
        "_detail": detail,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gt-root", required=True,
                    help="root dir holding per-sample ground_truth subdirs")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    sample_dirs = sorted(
        os.path.dirname(p) for p in glob.glob(
            os.path.join(a.gt_root, "*", "blastn_imgvr_raw.tsv")))
    if not sample_dirs:
        sys.exit(f"ERROR: no blastn_imgvr_raw.tsv found under {a.gt_root}")
    print(f"Auditing {len(sample_dirs)} sample(s) with E1b raw output...")

    rows, all_detail = [], []
    for sd in sample_dirs:
        r = audit_sample(sd)
        if r is None:
            continue
        all_detail.extend(r.pop("_detail"))
        rows.append(r)
        print(f"  {r['sample']:<14} E1b orig={r['n_e1b_orig']:>5} "
              f"corr={r['n_e1b_corrected']:>5} (changed {r['n_e1b_changed']:>4}: "
              f"+{r['n_e1b_gained']}/-{r['n_e1b_lost']})  "
              f"concord={r['frozen_concordance_pct']:.2f}%  "
              f"tier_changed={r['n_tier_changed']} (Tier1={r['n_tier1_changed']})")

    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    sum_cols = ["n_e1b_orig", "n_e1b_corrected", "n_e1b_gained", "n_e1b_lost",
                "n_e1b_changed", "n_tier_changed", "n_tier1_changed"]
    cols = ["sample"] + sum_cols[:5] + ["frozen_concordance_pct"] + sum_cols[5:]
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
        tot = {c: sum(r[c] for r in rows) for c in sum_cols}
        tot["sample"] = "TOTAL"
        # report the worst (min) per-sample concordance in the TOTAL row
        tot["frozen_concordance_pct"] = round(min(r["frozen_concordance_pct"] for r in rows), 3) if rows else 0.0
        w.writerow(tot)

    if all_detail:
        det_path = a.out.replace(".tsv", "_tier_changes.tsv")
        with open(det_path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(all_detail[0].keys()),
                               delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(all_detail)
        print(f"  wrote per-contig tier changes: {det_path}")

    tot_e1b = sum(r["n_e1b_changed"] for r in rows)
    tot_tier = sum(r["n_tier_changed"] for r in rows)
    tot_t1 = sum(r["n_tier1_changed"] for r in rows)
    print("\n" + "=" * 64)
    print("E1b PARSER AUDIT SUMMARY")
    print(f"  Total E1b call changes : {tot_e1b}")
    print(f"  Total TIER changes     : {tot_tier}")
    print(f"  Total Tier-1 changes   : {tot_t1}")
    print("-" * 64)
    if tot_tier == 0:
        print("  DECISION: n_tier_changed == 0 -> footnote only "
              "(E1b parser bug has no ground-truth impact; E1 = E1a OR E1b).")
    else:
        print("  DECISION: n_tier_changed > 0 -> add Methods/Supplement disclosure;"
              "\n            confirm Table S11 (E1b-excluded) bounds the effect."
              + ("\n            Tier-1 membership changes -> ESCALATE before submission."
                 if tot_t1 else ""))
    print("=" * 64)
    print(f"\nWrote {a.out}")


if __name__ == "__main__":
    main()
