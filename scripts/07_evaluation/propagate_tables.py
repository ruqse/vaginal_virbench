#!/usr/bin/env python3
"""Propagate the LOCKED 13-sample pooled multi-evidence + cascade numbers into the
results/tables files. Pure TSV reformatting (stdlib csv only) from:
  results/test_real/secondary_benchmark_pooled13/{overall_metrics.tsv,cascade/*}
Updates: Table 3 (multi-evidence cols), Table S12, Table S11, Table S13, Table S22.
The single-sample originals are preserved in git history.
"""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
POOL = ROOT / "results/test_real/secondary_benchmark_pooled13"
CAS = POOL / "cascade"
TAB = ROOT / "results/tables"

def read_tsv(p, comment="#"):
    with open(p) as fh:
        return [l.rstrip("\n") for l in fh]

# ---- load locked pooled metrics ----
pooled = {}
with open(POOL / "overall_metrics.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        pooled[r["tool"]] = r

# Table 3 is written in full by publish_parser_revision_tables.py (Track C and
# multi-evidence columns from their sources); it is not read back or patched here.

# ============================================================ Table S12
# tool, TP, FP, TN, FN, precision, recall, F1, MCC, AUPRC, n_predicted, n_total
s12 = TAB / "table_s12_secondary_benchmark_overall.tsv"
cols = ["tool", "TP", "FP", "TN", "FN", "precision", "recall", "F1", "MCC",
        "AUPRC", "n_predicted", "n_total"]
order = sorted(pooled.values(), key=lambda r: float(r["MCC"]), reverse=True)
with open(s12, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(cols)
    for r in order:
        w.writerow([r[c] for c in cols])
print(f"Table S12 updated ({len(order)} tools, n_total={order[0]['n_total']})")

# ============================================================ Table S11
# add deltas vs primary
s11_in = list(csv.DictReader(open(CAS / "s11_multi_evidence_sensitivity.tsv"),
                             delimiter="\t"))
out_cols = ["tool", "MCC_primary", "rank_primary",
            "MCC_no_E1b", "rank_no_E1b", "delta_MCC_no_E1b", "delta_rank_no_E1b",
            "MCC_merged_E1E2", "rank_merged_E1E2", "delta_MCC_merged_E1E2", "delta_rank_merged_E1E2",
            "MCC_no_E5", "rank_no_E5", "delta_MCC_no_E5", "delta_rank_no_E5"]
with open(TAB / "table_s11_multi_evidence_sensitivity.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(out_cols)
    for r in s11_in:
        mp, rp = float(r["MCC_primary"]), int(r["rank_primary"])
        row = [r["tool"], r["MCC_primary"], r["rank_primary"]]
        for v in ["no_E1b", "merged_E1E2", "no_E5"]:
            m, rk = float(r[f"MCC_{v}"]), int(r[f"rank_{v}"])
            row += [r[f"MCC_{v}"], r[f"rank_{v}"], f"{m-mp:.4f}", str(rk-rp)]
        w.writerow(row)
print(f"Table S11 updated (geNomad rank: "
      f"{[r['rank_primary'] for r in s11_in if r['tool']=='geNomad']})")

# ============================================================ Table S13
s13_in = list(csv.DictReader(open(CAS / "q2_prophage_detection.tsv"), delimiter="\t"))
FAMILY = {"geNomad": "DNN + markers", "VirSorter2": "Feature-based ML",
          "VirSorter": "Feature-based ML", "VIBRANT": "Feature-based ML",
          "VirFinder": "Feature-based ML", "DeepVirFinder": "Deep learning",
          "PPR-Meta": "Deep learning", "HVSeeker": "Deep learning",
          "Jaeger": "Deep learning", "Seeker": "Deep learning",
          "ViraLM": "Attention/transformer", "TransGINmer": "Attention/transformer",
          "MetaPhinder": "Reference-based", "Sourmash": "Reference-based"}
with open(TAB / "table_s13_q2_prophage_detection.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["tool", "type", "method_approach", "prophage_tp", "prophage_total",
                "prophage_recall", "freephage_pos", "freephage_total", "freephage_rate",
                "included_in_group_mean"])
    # Q4 group means exclude Jaeger/Seeker (contig-length incompatibility) and
    # Sourmash (degenerate: no positive call on either contig class).
    Q4_EXCLUDE = {"Jaeger", "Seeker", "Sourmash"}
    for r in s13_in:
        w.writerow([r["tool"], r["type"], FAMILY.get(r["tool"], ""),
                    r["prophage_tp"], r["prophage_total"], r["prophage_recall"],
                    r["freephage_pos"], r["freephage_total"], r["freephage_rate"],
                    r.get("included_in_group_mean",
                          "no" if r["tool"] in Q4_EXCLUDE else "yes")])
print(f"Table S13 updated ({len(s13_in)} tools; prophage_total=49)")

# ============================================================ Table S22
panelA = list(csv.DictReader(open(CAS / "cross_track_panelA.tsv"), delimiter="\t"))
def load_panel(name):
    return list(csv.DictReader(open(CAS / name), delimiter="\t"))
B = load_panel("cross_track_panelB.tsv")
C = load_panel("cross_track_panelC.tsv")
contexts = ["Track A (L1500 MCC)", "Track A (L3000 MCC)",
            "Track B (10x coverage, CST-I MCC)",
            "Multi-evidence benchmark (MCC)", "Track C (MCC)"]
with open(TAB / "table_s22_cross_track_spearman.tsv", "w", newline="") as fh:
    fh.write("# Supplementary Table S22. Cross-track Spearman rank correlations.\n")
    fh.write("# Panel A: per-tool MCC at each context (n = 14 tools). The multi-evidence\n")
    fh.write("#   column is the 13-sample pooled benchmark (>=1500 bp; Gate 0).\n")
    fh.write("# Panel B: pairwise Spearman (rho) with BH-adjusted p across 10 comparisons.\n")
    fh.write("# Panel C: restricted to the five top Track A L1500 tools.\n")
    fh.write("# Track B uses 10x CST-I (UC115_V2) MCC.\n\n")
    fh.write("## Panel A: per-tool MCC by context\n")
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["tool"] + contexts)
    for r in panelA:
        w.writerow([r["tool"]] + [r[c] for c in contexts])
    fh.write("\n## Panel B: pairwise Spearman rank correlations\n")
    w.writerow(["context_1", "context_2", "n_tools", "spearman_rho",
                "spearman_p", "spearman_p_BH"])
    for r in B:
        w.writerow([r["context_1"], r["context_2"], r["n_tools"],
                    r["spearman_rho"], r["spearman_p"], r["spearman_p_BH"]])
    fh.write("\n## Panel C: restricted to top-5 Track A L1500 tools "
             "(geNomad, ViraLM, VirSorter2, TransGINmer, VirSorter)\n")
    w.writerow(["context_1", "context_2", "n_tools", "spearman_rho",
                "spearman_p", "spearman_p_BH"])
    for r in C:
        w.writerow([r["context_1"], r["context_2"], r["n_tools"],
                    r["spearman_rho"], r["spearman_p"], r["spearman_p_BH"]])
print("Table S22 updated (Panel A/B/C)")
print("\nAll tables propagated.")
