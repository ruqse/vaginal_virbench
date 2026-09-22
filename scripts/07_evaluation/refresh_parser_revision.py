#!/usr/bin/env python3
"""Reproduce the September 2026 scoring revision from existing tool outputs.

Run independent stages with the SciPy-bundle/2024.05-gfbf-2024a module:
tracka, primary, trackb, trackc, rca. Then run mitch, bootstrap, publish.
No sequence classifier or assembly is rerun. Command logs and status are saved.
"""
from __future__ import annotations
import argparse
import csv
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
REV = ROOT / "results/regeneration_logs"  # command logs + stage status (gitignored)
TAB = ROOT / "results/tables"
UC = "UC028_V2 UC055_V1 UC055_V2 UC062_V2 UC065_V2 UC074_V2 UC084_V2 UC093_V2 UC093_V3 UC096_V2 UC115_V2 UC139_V2 UC164_V2".split()
RCA_ORDER = "UC084_V2 UC115_V2 UC164_V2 UC093_V2 UC139_V2 UC055_V1 UC065_V2 UC096_V2 UC028_V2 UC074_V2 UC093_V3 UC055_V2 UC062_V2".split()
stage = None


def run(script, *args):
    cmd = [sys.executable, script, *map(str, args)]
    log = REV / "logs" / f"{stage}_{Path(script).stem}_{int(time.time_ns())}.log"
    log.parent.mkdir(parents=True, exist_ok=True)
    print(json.dumps({"running": cmd, "log": str(log.relative_to(ROOT))}), flush=True)
    start = time.time()
    with log.open("w") as fh:
        result = subprocess.run(cmd, cwd=ROOT, stdout=fh, stderr=subprocess.STDOUT)
    with (REV / f"commands_{stage}.jsonl").open("a") as fh:
        fh.write(json.dumps({"command": cmd, "seconds": time.time() - start,
                             "returncode": result.returncode,
                             "log": str(log.relative_to(ROOT))}) + "\n")
    if result.returncode:
        print(log.read_text()[-8000:])
        raise RuntimeError(f"Failed: {cmd}; see {log}")


def read(path):
    with Path(path).open() as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))


def write(path, rows, fields=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields or list(rows[0]), delimiter="\t",
                           extrasaction="ignore", lineterminator="\n")
        w.writeheader(); w.writerows(rows)


def tracka():
    base = ROOT / "results/spike_in_benchmark"
    manifest = "data/spike_in/fragments/fragment_manifest.tsv"
    combined = []
    for length in [500, 1000, 1500, 3000, 5000, 10000]:
        d = base / f"L{length}_fragments"
        run("scripts/07_evaluation/benchmark_fragments.py", "--manifest", manifest,
            "--results-dir", d, "--output", d / "benchmark_results.tsv")
        rows = read(d / "benchmark_results.tsv")
        combined += [{"length_bin": f"L{length}", "length_bp": length, **r} for r in rows]
        if length == 1500:
            write(TAB / "table_s34_track_a_overall_metrics_full.tsv", sorted(rows, key=lambda r: -float(r["MCC"])))
            write(TAB / "table2_track_a_overall_metrics_L1500.tsv", rows,
                  ["tool", "precision", "recall", "F1", "MCC", "AUPRC"])
            shutil.copy2(d / "benchmark_results_by_category.tsv", TAB / "table_s1_track_a_by_category_L1500.tsv")
            shutil.copy2(d / "benchmark_results_by_source_genome.tsv", TAB / "table_s20_track_a_by_source_genome_L1500.tsv")
        else:
            write(ROOT / f"results/secondary_tables/table_s_track_a_L{length}.tsv", rows)
    write(TAB / "table_s7_track_a_overall_metrics_by_length.tsv", combined)
    run("scripts/07_evaluation/dump_pr_curve_scores.py", "--manifest", manifest,
        "--results-root", base, "--lengths", "L1500", "L3000",
        "--output", base / "pr_curve_scores.tsv.gz")
    novel = ROOT / "data/novel_spike_discovery/full_scale"
    for length in [1500, 3000]:
        run("scripts/07_evaluation/benchmark_fragments.py", "--manifest", novel / "fragments/fragment_manifest.tsv",
            "--results-dir", novel / f"tool_results/L{length}", "--output", novel / f"novel_benchmark_L{length}.tsv")
    run("scripts/07_evaluation/stratify_novel_panel.py")
    run("scripts/tables/build_table_s37.py")


def calibration():
    run("scripts/07_evaluation/calibrate_tool_thresholds.py",
        "--manifest", "data/spike_in/fragments/fragment_manifest.tsv",
        "--results-root", "results/spike_in_benchmark", "--lengths", "L1500_fragments", "L3000_fragments",
        "--output", "results/spike_in_benchmark/calibration/threshold_calibration.tsv")
    src = ROOT / "results/spike_in_benchmark/calibration/threshold_calibration.tsv"
    shutil.copy2(src, TAB / "table_s41_threshold_calibration.tsv")
    (ROOT / "results/secondary_tables").mkdir(parents=True, exist_ok=True); shutil.copy2(src, ROOT / "results/secondary_tables/table_s20_threshold_calibration.tsv")


def tracka_stats():
    run("scripts/07_evaluation/statistical_tests.py", "--manifest", "data/spike_in/fragments/fragment_manifest.tsv",
        "--results-root", "results/spike_in_benchmark", "--output-dir", "results/spike_in_benchmark/statistical_tests")
    for source, target in [("bootstrap_ci.tsv", "table_s4_bootstrap_confidence_intervals.tsv"),
                           ("cochrans_q.tsv", "table_s5_cochrans_q_test.tsv"),
                           ("mcnemar_significant.tsv", "table_s6_mcnemar_significant_pairs.tsv")]:
        shutil.copy2(ROOT / "results/spike_in_benchmark/statistical_tests" / source, TAB / target)


def primary():
    for suffix, gt in [("", "ground_truth_with_kraken2.tsv"), ("_no_e1b", "ground_truth_no_e1b.tsv"),
                       ("_merged_e1e2", "ground_truth_merged_e1e2.tsv"), ("_no_e5", "ground_truth_no_e5.tsv")]:
        run("scripts/07_evaluation/pooled_multievidence_benchmark.py", "--samples", *UC,
            "--results-dir", "results/test_real/full_run", "--gt-root", "results/test_real/ground_truth",
            "--gt-name", gt, "--min-length", "1500",
            "--output-dir", f"results/test_real/secondary_benchmark_pooled13{suffix}")


def trackb():
    run("scripts/07_evaluation/benchmark_track_b.py", "--assemblies-dir", "data/spike_in/assemblies",
        "--results-dir", "results/spike_in_benchmark", "--output-dir", "results/spike_in_benchmark/trackB_metrics")
    for source, target in [("trackB_metrics.tsv", "table_s8_track_b_metrics.tsv"),
                           ("trackB_recall_by_coverage.tsv", "table_s10_track_b_recall_by_coverage.tsv")]:
        shutil.copy2(ROOT / "results/spike_in_benchmark/trackB_metrics" / source, TAB / target)
    run("scripts/09_expansion/08_score_and_aggregate_track_b.py")
    run("scripts/09_expansion/11b_aggregate_excl.py")


def trackc():
    base = "results/test_real/coassembly"
    run("scripts/03_assembly/evaluate_rca_benchmark.py", "--ground-truth", f"{base}/rca_ground_truth.tsv",
        "--results-dir", "results/test_real/full_run/RCA_master", "--output-dir", f"{base}/benchmark")
    common = ["--coassembly-dir", base, "--checkv-dir", f"{base}/checkv",
              "--kraken2-output", f"{base}/kraken2/kraken2_output.tsv", "--kraken2-report", f"{base}/kraken2/kraken2_report.txt",
              "--results-dir", "results/test_real/full_run/RCA_master"]
    run("scripts/07_evaluation/track_c_threshold_sensitivity.py", *common, "--output-dir", f"{base}/sensitivity_R")
    run("scripts/07_evaluation/track_c_gmm_sensitivity.py", *common, "--rca-ground-truth", f"{base}/rca_ground_truth.tsv",
        "--output-dir", f"{base}/sensitivity_R_gmm")


def rca():
    run("scripts/07_evaluation/rca_concordance.py", "--samples", *RCA_ORDER,
        "--evidence-root", "results/test_real/ground_truth", "--results-root", "results/test_real/full_run",
        "--ground-truth-root", "results/test_real/ground_truth", "--e1a-root", "results/test_real/ground_truth",
        "--output-dir", "results/test_real/rca_concordance")
    run("scripts/07_evaluation/rca_sensitivity_grid.py", "--samples", *RCA_ORDER,
        "--stats-root", "results/test_real/ground_truth", "--results-root", "results/test_real/full_run",
        "--output-dir", "results/test_real/rca_sensitivity_grid", "--data-driven")


def mitch():
    run("scripts/07_evaluation/pooled_multievidence_benchmark.py", "--samples", *[f"MITCH{i:02d}" for i in range(1,31)],
        "--results-dir", "results/test_real/full_run", "--gt-root", "results/test_real/ground_truth",
        "--gt-name", "ground_truth_4line_kraken2.tsv", "--min-length", "1500", "--min-tier", "2",
        "--output-dir", "results/test_real/viralm_investigation/mitch_scaleup30")
    run("scripts/10_vmgc_crosscheck/scaleup_analysis.py")
    run("scripts/10_vmgc_crosscheck/scaleup_diversity_contrast.py")
    run("scripts/10_vmgc_crosscheck/06_export_mitch_supp_tables.py")
    run("scripts/07_evaluation/replication_sensitivity.py")
    run("scripts/10_vmgc_crosscheck/mitch_prevalence_diversity.py")


def bootstrap():
    run("scripts/07_evaluation/bootstrap_ci_trackc_multievidence.py")
    run("scripts/07_evaluation/revision_pairwise_bootstrap.py")


def publish():
    run("scripts/07_evaluation/recompute_cascade.py")
    run("scripts/07_evaluation/propagate_tables.py")
    run("scripts/07_evaluation/publish_parser_revision_tables.py")
    run("scripts/tables/build_table_s38.py")
    run("scripts/tables/build_main_tables.py")
    run("scripts/tables/build_supplementary_workbook.py")


if __name__ == "__main__":
    choices = ["tracka", "tracka_stats", "calibration", "primary", "trackb", "trackc", "rca", "mitch", "bootstrap", "publish"]
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("stage", choices=choices)
    stage = ap.parse_args().stage
    REV.mkdir(parents=True, exist_ok=True)
    (REV / f"status_{stage}.json").write_text(json.dumps({"status": "running", "started": time.time()}) + "\n")
    try:
        globals()[stage]()
    except BaseException as exc:
        (REV / f"status_{stage}.json").write_text(json.dumps({"status": "failed", "failed": time.time(), "error": str(exc)}) + "\n")
        raise
    (REV / f"status_{stage}.json").write_text(json.dumps({"status": "complete", "completed": time.time()}) + "\n")
