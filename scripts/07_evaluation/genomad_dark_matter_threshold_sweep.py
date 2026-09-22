#!/usr/bin/env python3
"""
genomad_dark_matter_threshold_sweep.py

Re-score geNomad's virus_score output for the 48 Track C dark-matter contigs at
a grid of thresholds. We ran geNomad in the main benchmark at its default
threshold (virus_score >= 0.7; THRESHOLDS.md); all other tools ran at 0.5.
That threshold asymmetry could inflate geNomad's low dark-matter recall
(1/48, 2%). This script re-counts detection at 0.3 / 0.5 / 0.6 / 0.7 / 0.8 /
0.9 to separate "genuinely blind on novel short contigs" from "conservative
threshold penalises short novel contigs."

Inputs (pre-existing):
  - results/.../RCA_master/genomad/.../RCA_master_contigs_aggregated_classification.tsv
  - results/test_real/coassembly/rca_dark_matter_contig_ids.txt
  - results/tables/table_s14_dark_matter_per_contig.tsv  (for length + R per contig)

Output:
  - results/tables/table_s_genomad_dark_matter_threshold_sweep.tsv
  - stdout summary
"""

import csv
from pathlib import Path

PROJ = Path(__file__).resolve().parents[2]
GENOMAD_TSV = PROJ / ("results/test_real/full_run/RCA_master/genomad/"
                     "RCA_master_contigs_aggregated_classification/"
                     "RCA_master_contigs_aggregated_classification.tsv")
DARK_IDS = PROJ / "results/test_real/coassembly/rca_dark_matter_contig_ids.txt"
PER_CONTIG = PROJ / "results/tables/table_s14_dark_matter_per_contig.tsv"
OUT_TSV = PROJ / "results/tables/table_s17_genomad_dark_matter_threshold_sweep.tsv"

THRESHOLDS = [0.30, 0.50, 0.60, 0.70, 0.80, 0.90]
DEFAULT_THRESHOLD = 0.70


def main():
    dark_ids = [l.strip() for l in DARK_IDS.read_text().splitlines() if l.strip()]
    print(f"dark-matter contigs: {len(dark_ids)}")

    # Load per-contig metadata (length, R, current detection count) from Table S11
    per_contig = {}
    with open(PER_CONTIG) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            per_contig[row["contig_id"]] = row

    # Load geNomad virus scores for just the dark-matter contigs
    scores = {}
    with open(GENOMAD_TSV) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            cid = row["seq_name"]
            if cid in set(dark_ids):
                scores[cid] = float(row["virus_score"])
    print(f"geNomad scores found for: {len(scores)}/{len(dark_ids)}")

    # Summarise score distribution on the 48
    sorted_s = sorted(scores.values(), reverse=True)
    print(f"\nScore distribution on 48 dark-matter contigs:")
    print(f"  max = {sorted_s[0]:.4f}")
    print(f"  top-5: {[f'{s:.4f}' for s in sorted_s[:5]]}")
    print(f"  median = {sorted_s[len(sorted_s)//2]:.4f}")
    print(f"  min = {sorted_s[-1]:.4f}")

    # Threshold sweep
    rows = []
    print(f"\nDetection at varying thresholds:")
    for t in THRESHOLDS:
        detected = [c for c, s in scores.items() if s >= t]
        n = len(detected)
        is_default = (t == DEFAULT_THRESHOLD)
        marker = " (DEFAULT)" if is_default else ""
        print(f"  virus_score >= {t:.2f}: {n:2d}/48 ({100*n/48:4.1f}%){marker}")
        rows.append({
            "virus_score_threshold": f"{t:.2f}",
            "n_detected": n,
            "detection_rate": f"{n/48:.4f}",
            "is_default": int(is_default),
            "detected_contig_ids": ";".join(sorted(detected)) if detected else "",
        })

    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_TSV, "w", newline="") as f:
        fields = ["virus_score_threshold", "n_detected", "detection_rate",
                  "is_default", "detected_contig_ids"]
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {OUT_TSV}")

    # Also list the per-contig scores for the 48 for appendix
    appendix_path = PROJ / "results/secondary_tables/table_s24_genomad_dark_matter_per_contig_scores.tsv"; appendix_path.parent.mkdir(parents=True, exist_ok=True)
    with open(appendix_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["contig_id", "contig_length", "enrichment_R",
                    "genomad_virus_score",
                    "detected_at_0.50", "detected_at_0.70"])
        for cid in dark_ids:
            s = scores.get(cid, float("nan"))
            meta = per_contig.get(cid, {})
            w.writerow([
                cid,
                meta.get("contig_length", ""),
                meta.get("enrichment_R", ""),
                f"{s:.4f}" if s == s else "NA",
                int(s >= 0.50) if s == s else "NA",
                int(s >= 0.70) if s == s else "NA",
            ])
    print(f"wrote {appendix_path}")


if __name__ == "__main__":
    main()
