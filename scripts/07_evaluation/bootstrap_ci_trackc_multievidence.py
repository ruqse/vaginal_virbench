#!/usr/bin/env python3
"""Bootstrap 95% CIs for the multi-evidence benchmark and Track C.

Uses the SAME resampling routine as Track A's Table S4
(statistical_tests.bootstrap_metric: nonparametric bootstrap over contigs,
1,000 resamples, percentile interval, seed 42), so intervals are comparable
across tracks.

Outputs:
  results/ci/multievidence_bootstrap_ci.tsv
  results/ci/track_c_bootstrap_ci.tsv
"""
import gzip
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "07_evaluation"))
sys.path.insert(0, str(ROOT / "scripts" / "03_assembly"))

from statistical_tests import bootstrap_metric  # noqa: E402

METRICS = ["MCC", "AUPRC", "precision", "recall", "F1"]
N_BOOT = 1000
SEED = 42


def rows_for(tool, y_true, y_pred, y_scores, extra):
    out = []
    for m in METRICS:
        b = bootstrap_metric(y_true, y_pred, y_scores, m, N_BOOT, SEED, 0.05)
        out.append({"tool": tool, "metric": m, **extra, **b})
    return out


def multievidence():
    f = ROOT / "results/test_real/secondary_benchmark_pooled13/pooled_tool_predictions.tsv.gz"
    per_tool = {}
    with gzip.open(f, "rt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        i = {c: n for n, c in enumerate(hdr)}
        for line in fh:
            p = line.rstrip("\n").split("\t")
            t = p[i["tool"]]
            d = per_tool.setdefault(t, {"y": [], "p": [], "s": []})
            d["y"].append(1 if p[i["label"]] == "viral" else 0)
            d["p"].append(int(float(p[i["y_pred"]])))
            sc = p[i["score"]]
            d["s"].append(float(sc) if sc not in ("", "NA", "nan") else 0.0)

    rows = []
    for tool, d in sorted(per_tool.items()):
        y_true = np.array(d["y"]); y_pred = np.array(d["p"]); y_sc = np.array(d["s"])
        rows += rows_for(tool, y_true, y_pred, y_sc,
                         {"n_total": len(y_true), "n_pos": int(y_true.sum())})
        print(f"  multievidence {tool:14s} n={len(y_true)} pos={int(y_true.sum())}")
    return rows


def track_c():
    import evaluate_rca_benchmark as ev
    gt = ev.load_rca_ground_truth(str(ROOT / "results/test_real/coassembly/rca_ground_truth.tsv"))
    gt_ids = set(gt)
    scores = ev.discover_tool_outputs(
        Path(ROOT / "results/test_real/full_run/RCA_master"), "RCA_master", gt_ids=gt_ids)

    contig_ids = sorted(gt)
    y_true = np.array([1 if gt[c]["label"] == "viral" else 0 for c in contig_ids])

    rows = []
    for tool in sorted(scores):
        ts = scores[tool]
        thr = ev.TOOL_THRESHOLDS.get(tool, 0.5)
        pthr = ev.TOOL_PVALUE_THRESHOLDS.get(tool)
        sample_val = next((v for v in ts.values()), None)
        if isinstance(sample_val, tuple):
            y_sc = np.array([ts[c][0] if c in ts else 0.0 for c in contig_ids])
            y_pv = np.array([ts[c][1] if c in ts else 1.0 for c in contig_ids])
            y_pred = (((y_sc >= thr) & (y_pv < pthr)) if pthr is not None
                      else (y_sc >= thr)).astype(int)
        else:
            y_sc = np.array([float(ts.get(c, 0.0)) for c in contig_ids])
            y_pred = (y_sc >= thr).astype(int)
        rows += rows_for(tool, y_true, y_pred, y_sc,
                         {"n_total": len(y_true), "n_pos": int(y_true.sum())})
        print(f"  track_c {tool:14s} n={len(y_true)} pos={int(y_true.sum())}")
    return rows


def write(rows, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    cols = ["tool", "metric", "n_total", "n_pos", "point", "ci_lo", "ci_hi", "se"]
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "NA")) for c in cols) + "\n")
    print(f"wrote {path} ({len(rows)} rows)")


if __name__ == "__main__":
    print("== multi-evidence ==")
    write(multievidence(), ROOT / "results/ci/multievidence_bootstrap_ci.tsv")
    print("== track C ==")
    write(track_c(), ROOT / "results/ci/track_c_bootstrap_ci.tsv")
