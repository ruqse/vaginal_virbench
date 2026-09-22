#!/usr/bin/env python3
"""Publish regenerated numeric outputs while preserving established table schemas."""
import csv
import json
import re
import shutil
from pathlib import Path
from refresh_parser_revision import ROOT, TAB, read, write
from metric_display import format_mcc_from_counts

TEMPLATES = TAB / "templates"

# propagate_tables refreshes Panels A-C; the independent-cohort panel has its
# own reproducible analysis and must survive rebuilding the combined S22.
s22 = TAB / "table_s22_cross_track_spearman.tsv"
body = s22.read_text().split("## Panel D:", 1)[0].rstrip()
s22.write_text(body + "\n\n## Panel D: independent-cohort rank replication and level agreement\n" +
               (TAB / "table_s22d_replication_sensitivity.tsv").read_text())
paired = json.loads((ROOT / "results/statistics/paired_mcc_intervals.json").read_text())
write(TAB / "table_s27_pairwise_mcc.tsv", [{"context": context, **row}
      for context, rows in paired.items() for row in rows])


def attach_ci(rows, path):
    intervals = {(r["tool"], r["metric"]): r for r in read(path)}
    for r in rows:
        for metric in ["MCC", "AUPRC"]:
            ci = intervals[r["tool"], metric]
            if r[metric] not in {"NA", "", "nan"}:
                assert abs(float(r[metric]) - float(ci["point"])) < .00011, (r["tool"], metric, r[metric], ci)
            r[f"{metric}_CI_low"] = ci["ci_lo"]
            r[f"{metric}_CI_high"] = ci["ci_hi"]
    return rows


me = read(ROOT / "results/test_real/secondary_benchmark_pooled13/overall_metrics.tsv")
metadata = {r["tool"]: r for r in read(TAB / "table_s21_tool_panel_metadata.tsv")}
me = attach_ci(me, ROOT / "results/ci/multievidence_bootstrap_ci.tsv")
fields = ["tool", "TP", "FP", "TN", "FN", "precision", "recall", "F1", "MCC", "MCC_CI_low", "MCC_CI_high",
          "AUPRC", "AUPRC_CI_low", "AUPRC_CI_high", "n_predicted", "n_total"]
write(TAB / "table_s12_secondary_benchmark_overall.tsv", sorted(me, key=lambda r: -float(r["MCC"])), fields)

tc_dir = ROOT / "results/test_real/coassembly"
tc = read(tc_dir / "benchmark/overall_metrics.tsv")
tc = attach_ci(tc, ROOT / "results/ci/track_c_bootstrap_ci.tsv")
for r in tc:
    r["method_approach"] = metadata[r["tool"]]["method_approach"]
    r["scope"] = metadata[r["tool"]]["scope_label"]
fields = ["tool", "method_approach", "scope", "TP", "FP", "TN", "FN", "precision", "recall", "FDR", "F1", "MCC",
          "MCC_CI_low", "MCC_CI_high", "AUPRC", "AUPRC_CI_low", "AUPRC_CI_high", "n_predicted", "n_total"]
write(TAB / "table_s35_track_c_overall_metrics_full.tsv", sorted(tc, key=lambda r: -float(r["MCC"])), fields)
me_by = {r["tool"]: r for r in me}
write(TAB / "table3_track_c_overall_metrics.tsv", [{"tool": r["tool"], "track_c_mcc": r["MCC"],
    "track_c_auprc": r["AUPRC"], "multi_evidence_mcc": me_by[r["tool"]]["MCC"],
    "multi_evidence_auprc": me_by[r["tool"]]["AUPRC"]} for r in sorted(tc, key=lambda r: -float(r["MCC"]))])

copies = {
    "benchmark/dark_matter_per_contig.tsv": "table_s14_dark_matter_per_contig.tsv",
    "benchmark/dark_matter_detection_rates.tsv": "table_s15_dark_matter_detection_rates.tsv",
    "sensitivity_R/track_c_sensitivity_mcc.tsv": "table_s18_track_c_R_sensitivity_mcc.tsv",
    "sensitivity_R/track_c_sensitivity_ranking.tsv": "table_s18_track_c_R_sensitivity_ranking.tsv",
}
for source, dest in copies.items():
    shutil.copy2(tc_dir / source, TAB / dest)

fit = read(tc_dir / "sensitivity_R_gmm/track_c_gmm_fit.tsv")[0]
counts = read(tc_dir / "sensitivity_R_gmm/track_c_gmm_counts.tsv")[0]
gmm = read(tc_dir / "sensitivity_R_gmm/track_c_gmm_vs_default.tsv")
out = TAB / "table_s23_track_c_gmm_sensitivity.tsv"
write(out, gmm)
body = out.read_text()
header = ["# Table S23. Track C performance under Gaussian-mixture relabelling.",
          "# Final reported viral calls are mapped to known parent contigs.",
          "# GMM fitting inputs and label rules are unchanged by the parser correction."]
header += [f"# fit {key} = {value}" for key, value in fit.items()]
header += [f"# counts {key} = {value}" for key, value in counts.items()]
out.write_text("\n".join(header) + "\n" + body)

for source, dest, keys in [
    ("concordance_by_cst.tsv", "table_s2_track_c_by_cst.tsv", ["cst", "tool"]),
    ("concordance_by_virus_category.tsv", "table_s3_track_c_by_virus_category.tsv", ["sample", "tool", "virus_category"]),
]:
    # Row order, keys and author-owned columns come from a template whose numeric
    # cells are blank; the published table is never read back, so a stale number
    # cannot survive a rebuild.
    old = read(TEMPLATES / dest.replace(".tsv", ".template.tsv"))
    updated = {tuple(r[k] for k in keys): r for r in read(ROOT / "results/test_real/rca_concordance" / source)}
    assert {tuple(r[k] for k in keys) for r in old} == set(updated), dest
    for r in old:
        new = updated[tuple(r[k] for k in keys)]
        filled = set(keys)
        for k in r:
            if k in new:
                r[k] = new[k]
                filled.add(k)
        if "method_approach" in r:
            r["method_approach"] = metadata[r["tool"]]["method_approach"]
            filled.add("method_approach")
        if "tool_scope" in r:
            r["tool_scope"] = metadata[r["tool"]]["scope_label"]
            filled.add("tool_scope")
        unfilled = [k for k, v in r.items() if v == "" and k not in filled]
        assert not unfilled, (dest, keys, unfilled)
    write(TAB / dest, old)

for source, dest in [("grid_results.tsv", "table_s19_rca_sensitivity_grid.tsv"),
                     ("grid_ranking_stability.tsv", "table_s19_rca_ranking_stability.tsv")]:
    shutil.copy2(ROOT / "results/test_real/rca_sensitivity_grid" / source, TAB / dest)

# The full recommendation matrix retains author-owned prose (from its template) but
# never stale MCCs: the template's MCC cells are blank and are recomputed here.
full = read(TEMPLATES / "table2_tool_recommendations_full.template.tsv")
ta_by = {r["tool"]: r for r in read(TAB / "table_s34_track_a_overall_metrics_full.tsv")}
tc_by = {r["tool"]: r for r in tc}
pro = {r["tool"]: r for r in read(TAB / "table_s13_q2_prophage_detection.tsv")}
for r in full:
    t = r["Tool"]
    if t != "Jaeger":
        r["Track A MCC (L1500)"] = format_mcc_from_counts(ta_by[t])
    r["Track C MCC (>=2500 bp)"] = format_mcc_from_counts(tc_by[t])
    if t == "VIBRANT":
        pct = 100 * int(pro[t]["prophage_tp"]) / int(pro[t]["prophage_total"])
        r["Tier"] = "Recommended"
        r["Research Goal"] = "High-precision phage identification on real assemblies"
        r["Min Contig Length"] = "1,000 bp; ≥4 predicted ORFs"
        r["Key Strength"] = (f"Highest MiTCH MCC; primary-assembly MCC difference from geNomad unresolved; "
                             f"{pct:.1f}% prophage recall; metabolic-gene annotation")
        r["Key Caveat"] = (f"Primary-assembly precision {float(me_by[t]['precision']):.3f} with recall {float(me_by[t]['recall']):.3f}; "
                            "targets prokaryotic viruses; zero eukaryote-infecting virus recall at 1,500 bp")
    elif t == "geNomad":
        r["Research Goal"] = "Balanced identification across evaluated designs"
        r["Key Strength"] = "Highest Track A MCC at 1,500 bp; primary-assembly MCC difference from VIBRANT unresolved"
    elif t == "ViraLM":
        r["Key Caveat"] = ("Higher FDR than geNomad (~6x at 1500 bp); catalogue use requires independent false-positive assessment. "
                            "CheckV assesses completeness and flanking host contamination; it does not confirm viral identity or remove false-positive viral predictions "
                            "(Nayfach et al. 2021; MIUViG quality reporting: Roux et al. 2019). Combined workflow not benchmarked here.")
for r in full:
    assert r["Track C MCC (>=2500 bp)"] and r["Track A MCC (L1500)"], r["Tool"]
write(TAB / "table2_tool_recommendations_full.tsv", sorted(full, key=lambda row: {"Recommended": 0, "Complementary": 1}.get(row["Tier"], 2)))
print("Published corrected S2/S3/S12/S14/S15/S18/S19/S23/S35 and source main-table metrics.")
