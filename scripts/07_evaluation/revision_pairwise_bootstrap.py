#!/usr/bin/env python3
"""Refresh paired MCC intervals using documented contig/sample bootstrap units."""
import csv
import gzip
import json
import sys
from itertools import combinations
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/03_assembly"))
import evaluate_rca_benchmark as ev
TOOLS = ["geNomad", "ViraLM", "VirSorter2", "VIBRANT", "Jaeger"]


def mcc(y, p):
    tp = int(np.sum((y == 1) & (p == 1))); fp = int(np.sum((y == 0) & (p == 1)))
    tn = int(np.sum((y == 0) & (p == 0))); fn = int(np.sum((y == 1) & (p == 0)))
    den = np.sqrt(float((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return float((tp*tn-fp*fn)/den) if den else 0.


def pooled(path):
    data = {}
    with gzip.open(path, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["tool"] in TOOLS:
                d = data.setdefault(r["contig_id"], {"sample": r["sample"], "y": int(r["label"] == "viral")})
                d[r["tool"]] = int(r["y_pred"])
    ids = sorted(data)
    return np.array([data[c]["y"] for c in ids]), {t: np.array([data[c][t] for c in ids]) for t in TOOLS}, np.array([data[c]["sample"] for c in ids])


def summarise(y, predictions, sample_labels=None):
    rng = np.random.default_rng(42 if sample_labels is None else 12345)
    pairs = list(combinations(TOOLS, 2))
    point = {t: mcc(y, predictions[t]) for t in TOOLS}
    diffs = {p: [] for p in pairs}
    if sample_labels is not None:
        samples = sorted(set(sample_labels))
        blocks = {s: np.flatnonzero(sample_labels == s) for s in samples}
    for _ in range(2000):
        if sample_labels is None:
            draw = rng.integers(0, len(y), len(y))
        else:
            draw = np.concatenate([blocks[samples[i]] for i in rng.integers(0, len(samples), len(samples))])
        vals = {t: mcc(y[draw], predictions[t][draw]) for t in TOOLS}
        for a, b in pairs:
            diffs[a, b].append(vals[a] - vals[b])
    result = []
    for a, b in pairs:
        vals = np.asarray(diffs[a, b]); lo, hi = np.percentile(vals, [2.5,97.5])
        result.append({"tool_a": a, "tool_b": b, "delta_a_minus_b": point[a]-point[b],
                       "ci_low": float(lo), "ci_high": float(hi), "ci_excludes_zero": bool(lo > 0 or hi < 0),
                       "bootstrap_draws": 2000, "unit": "contig" if sample_labels is None else "sample",
                       "seed": 42 if sample_labels is None else 12345})
    return result


y, p, sample = pooled(ROOT / "results/test_real/secondary_benchmark_pooled13/pooled_tool_predictions.tsv.gz")
out = {"multi_evidence": summarise(y, p)}
gt = ev.load_rca_ground_truth(str(ROOT / "results/test_real/coassembly/rca_ground_truth.tsv"))
scores = ev.discover_tool_outputs(ROOT / "results/test_real/full_run/RCA_master", "RCA_master", gt_ids=set(gt))
ids = sorted(gt); y = np.array([gt[c]["label"] == "viral" for c in ids], dtype=int)
p = {}
for t in TOOLS:
    threshold = ev.TOOL_THRESHOLDS[t]
    pvalue = ev.TOOL_PVALUE_THRESHOLDS.get(t)
    def call(c):
        v = scores[t].get(c, (0.,1.) if pvalue is not None else 0.)
        if isinstance(v, tuple):
            return int(v[0] >= threshold and (pvalue is None or v[1] < pvalue))
        return int(v >= threshold)
    p[t] = np.array([call(c) for c in ids])
out["track_c"] = summarise(y, p)
y, p, sample = pooled(ROOT / "results/test_real/viralm_investigation/mitch_scaleup30/pooled_tool_predictions.tsv.gz")
out["mitch"] = summarise(y, p, sample)
dest = ROOT / "results/statistics/paired_mcc_intervals.json"
dest.write_text(json.dumps(out, indent=2) + "\n")
print(dest)
