#!/usr/bin/env python3
"""Export the MITCH external-validation results as published supplementary
tables S27-S32.

The source of truth is computed upstream by ``scaleup_analysis.py`` and
``scaleup_diversity_contrast.py`` (per-tool ranking, per-CST/per-mgCST best
tools, per-sample MCC, diversity contrast) and by the mgCST assignment
pipeline (per-sample CST). This step copies those outputs verbatim into the
``results/tables/table_s<N>_*.tsv`` published names so every Table S27-S32 cited
in the manuscript resolves to a backing file. Values are never transformed;
only column selection (S32, to drop internal sequencing identifiers) differs.

Mapping:
  S27  scaleup30_ranking_ci.tsv          per-tool pooled MCC + sample-block CIs
  S28  scaleup30_per_CST.tsv             top-3 tools per traditional CST
  S29  scaleup30_per_mgCST.tsv           top-3 tools per metagenomic CST
  S30  scaleup30_per_sample_MCC.tsv      per-sample MCC matrix (backs Fig. S16)
  S31  scaleup30_diversity_contrast.tsv  high- vs low-diversity delta MCC + CI
  S32  mgCSTs_with_CSTs_mitch30.csv      per-sample CST / mgCST assignments

Run:
    python3 scripts/10_vmgc_crosscheck/06_export_mitch_supp_tables.py
"""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCALE = ROOT / "results/test_real/viralm_investigation/mitch_scaleup30"
MGCST = ROOT / "results/mgcst_mitch30/mgCSTs_with_CSTs_mitch30.csv"
OUT = ROOT / "results/tables"

# S27-S31: verbatim copies of the scaleup30 outputs (no numeric transformation)
VERBATIM = {
    "scaleup30_ranking_ci.tsv":         "table_s27_mitch_ranking_ci.tsv",
    "scaleup30_per_CST.tsv":            "table_s28_mitch_per_cst.tsv",
    "scaleup30_per_mgCST.tsv":          "table_s29_mitch_per_mgcst.tsv",
    "scaleup30_per_sample_MCC.tsv":     "table_s30_mitch_per_sample_mcc.tsv",
    "scaleup30_diversity_contrast.tsv": "table_s31_mitch_diversity_contrast.tsv",
}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for src, dest in VERBATIM.items():
        text = (SCALE / src).read_text(encoding="utf-8")
        (OUT / dest).write_text(text, encoding="utf-8")
        print(f"  {src} -> {dest}")

    # S32: per-sample CST assignments; drop internal token/ssid identifiers
    keep = ["sampleID", "CST", "CST_subtype", "CST_score", "mgCST", "mgCST_score"]
    with open(MGCST, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    dest = OUT / "table_s32_mitch_per_sample_cst.tsv"
    with open(dest, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(keep)
        for row in rows:
            w.writerow([row[k] for k in keep])
    print(f"  {MGCST.name} -> {dest.name}  ({len(rows)} samples)")


if __name__ == "__main__":
    main()
