#!/usr/bin/env python3
"""Regenerate Table S4 (Track A per-tool bootstrap CIs) with the tie-aware AUPRC.

Uses the identical code path, seed and parameters as statistical_tests.py so
the MCC/F1/precision/recall rows reproduce byte-for-byte; only the AUPRC rows
change, because _auprc no longer depends on input row order.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "07_evaluation"))

from statistical_tests import (build_prediction_arrays, bootstrap_metric,  # noqa: E402
                               TOOL_THRESHOLDS, parse_length_from_name)
from benchmark_fragments import load_manifest  # noqa: E402

MANIFEST = ROOT / "data/spike_in/fragments/fragment_manifest.tsv"
RESULTS = ROOT / "results/spike_in_benchmark"
OUT = ROOT / "results/tables/table_s4_bootstrap_confidence_intervals.tsv"
METRICS = ["MCC", "AUPRC", "F1", "precision", "recall"]

truth_all = load_manifest(str(MANIFEST))
dirs = sorted([d for d in RESULTS.iterdir()
               if d.is_dir() and parse_length_from_name(d.name)],
              key=lambda d: parse_length_from_name(d.name))

rows = []
for rd in dirs:
    L = parse_length_from_name(rd.name)
    truth = {f: i for f, i in truth_all.items() if i["fragment_length"] == L}
    pred = build_prediction_arrays(truth, rd, list(TOOL_THRESHOLDS.keys()), TOOL_THRESHOLDS)
    print(f"{rd.name}: {len(truth)} fragments, {len(pred)} tools", flush=True)
    for tool in TOOL_THRESHOLDS:
        if tool not in pred:
            continue
        y_true, y_pred, y_scores = pred[tool]
        for m in METRICS:
            b = bootstrap_metric(y_true, y_pred, y_scores, m, 1000, 42, 0.05)
            rows.append({"length_bin": rd.name, "length_bp": L, "tool": tool,
                         "metric": m, **b})

cols = ["length_bin", "length_bp", "tool", "metric", "point", "ci_lo", "ci_hi", "se"]
with open(OUT, "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for r in rows:
        fh.write("\t".join(str(r.get(c, "NA")) for c in cols) + "\n")
print(f"wrote {OUT} ({len(rows)} rows)")
